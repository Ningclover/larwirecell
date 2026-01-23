# Plan: Extend TrackIDPIDMap2h5 with Richer Per-TrackID Truth Info

## Context

The current `TrackIDPIDMap2h5.cxx` stores only: `track_ids`, `pids`, `mother_ids`.
Goal: add G4 process name, mother's PDG, energy fraction, and optionally cover ALL simulated
particles (not just those that ionized in TPC) to support classification of:
**Tracks, Showers, Blips, Michel electrons, Delta-rays** — as described in the slide.

---

## Q1: How are TrackIDs currently found?

**Source:** `sim::SimChannel` via `sc.TDCIDEMap()` → `sim::IDE::trackID`

- Only TrackIDs that **ionized in the TPC active volume** appear.
- TrackID = 0 is skipped.
- EM shower secondaries may have a modified trackID; `ide.origTrackID` tracks the generator-level ancestor.
- **No connection to reco_charge frames** — this component reads truth only, independent of reco imaging.

---

## Q2: Can we store ALL simulated TrackIDs (not just TPC ionizers)?

**Yes.** Two complementary approaches:

### Approach A: Read MCParticle directly from the art event

```cpp
art::Handle<std::vector<simb::MCParticle>> particleHandle;
event.getByLabel("largeant", particleHandle);
// This gives ALL stored Geant4 particles regardless of TPC ionization
```

MCParticle collection contains ALL G4-tracked particles **except** those dropped by:
- Energy threshold (`EnergyCut`, default 0.0 GeV — usually not an issue)
- Excluded physics processes when `keepEMShowerDaughters=false`

### Approach B: Also read droppedMCParticles (for completeness)

```cpp
art::Handle<std::vector<simb::MCParticle>> droppedHandle;
event.getByLabel("largeant", "droppedMCParticles", droppedHandle);
```

Also available: `sim::ParticleAncestryMap` — maps dropped TrackIDs → their first stored ancestor.

### Recommendation

For pixel-level labeling needs (Blips, Michel, Delta-rays), **Approach A** (full MCParticle list)
is sufficient and covers essentially all physics-relevant particles. Dropped particles are usually
very low-energy EM secondaries that don't produce detectable signals anyway.

---

## Q3: What new data should be stored?

### New HDF5 datasets (parallel arrays per frame_ident)

| Dataset | Type | Source | Purpose |
|---|---|---|---|
| `processes` | int32 | `particle->Process()` → enum | G4 process name as integer code |
| `mother_pids` | int32 | lookup mother's MCParticle → PdgCode | Parent particle's PDG |
| `energy_fracs` | float32 | sum `ide.energyFrac` across all IDEs | TPC ionization energy fraction |

### G4 process → integer mapping (following CellTree convention)

| Code | Process | Particle category |
|---|---|---|
| 0 | `primary` | Primary tracks/showers from neutrino vertex |
| 1 | `Decay` | Decay products — **Michel electrons** (if parent=muon) |
| 2 | `eIoni` | Electron ionization secondaries — **Delta-rays / Blips** |
| 3 | `muIoni` | Muon ionization secondaries — **Delta-rays** |
| 4 | `eBrem` | Bremsstrahlung gammas — shower component |
| 5 | `compt` | Compton scattering |
| 6 | `phot` | Photoelectric — **Blips** |
| 7 | `conv` | Pair production |
| 8 | `hIoni` | Hadron ionization |
| 9 | `nCapture` | Neutron capture |
| -1 | (unknown) | Fallback for unlisted processes |

### Classification logic (downstream, in Python/ML training):

| Label | Rule |
|---|---|
| **Track** | PDG ∈ {13, 211, 2212, 321, ...}, process = primary or Decay |
| **Shower** | PDG ∈ {11, 22}, process = primary, eBrem, compt, conv |
| **Michel** | PDG = 11 (e-), process = Decay, mother_pid = ±13 |
| **Delta-ray** | PDG = 11, process = eIoni or muIoni, mother_pid = muon/pion/proton |
| **Blip** | PDG = 11 or 22, process = phot or eIoni, low energy |

---

## Implementation Plan

### Step 1: Change TrackID source from SimChannel-only to MCParticle + SimChannel

In `visit()` (the art event hook in TrackIDPIDMap2h5.cxx):

```cpp
// EXISTING: iterate SimChannel for TPC ionizers
// ADD: also read all MCParticles
art::Handle<std::vector<simb::MCParticle>> particleHandle;
event.getByLabel(m_particle_label, particleHandle);  // "largeant"
for (auto const& particle : *particleHandle) {
    int tid = particle.TrackId();
    m_trackid_to_pid[tid]       = particle.PdgCode();
    m_trackid_to_motherid[tid]  = particle.Mother();
    m_trackid_to_process[tid]   = processCode(particle.Process());
    // mother_pid requires lookup below
}
// Then also iterate SimChannel to get energy_fracs
for (auto const& sc : *simchan_handle) {
    for (auto const& [tdc, ides] : sc.TDCIDEMap()) {
        for (auto const& ide : ides) {
            if (ide.trackID == 0) continue;
            m_trackid_to_energyfrac[ide.trackID] += ide.energyFrac;
        }
    }
}
```

### Step 2: Add mother_pid lookup

After building the trackid→pid map from MCParticle, fill mother_pid:
```cpp
for (auto& [tid, motherid] : m_trackid_to_motherid) {
    if (motherid != 0 && m_trackid_to_pid.count(motherid)) {
        m_trackid_to_motherpid[tid] = m_trackid_to_pid[motherid];
    }
}
```

### Step 3: Add new member variables

In `TrackIDPIDMap2h5.h`:
```cpp
std::unordered_map<int, int>   m_trackid_to_process;
std::unordered_map<int, int>   m_trackid_to_motherpid;
std::unordered_map<int, float> m_trackid_to_energyfrac;
std::string m_particle_label{"largeant"};  // FCL configurable
```

### Step 4: Write new HDF5 datasets

In the existing write loop (where `track_ids`, `pids`, `mother_ids` are written):
```cpp
// After sorting track_ids vector:
std::vector<int>   processes(n), mother_pids(n);
std::vector<float> energy_fracs(n);
for (size_t i = 0; i < n; ++i) {
    int tid = track_ids[i];
    processes[i]    = m_trackid_to_process.count(tid) ? m_trackid_to_process[tid] : -1;
    mother_pids[i]  = m_trackid_to_motherpid.count(tid) ? m_trackid_to_motherpid[tid] : 0;
    energy_fracs[i] = m_trackid_to_energyfrac.count(tid) ? m_trackid_to_energyfrac[tid] : 0.0f;
}
// Write to HDF5 group alongside existing datasets
```

---

## Files to Modify

| File | Change |
|---|---|
| `TrackIDPIDMap2h5.h` | Add new member maps, `m_particle_label` config |
| `TrackIDPIDMap2h5.cxx` | Add MCParticle reading, new maps, new HDF5 writes |

## Files for Reference

| File | Purpose |
|---|---|
| `larreco/WireCell/CellTree_module.cc` | Process map enum, MCParticle access pattern |
| `lardataobj/Simulation/SimChannel.h` | IDE struct: `energyFrac`, `origTrackID` fields |
| `larg4/pluginActions/ParticleListAction.cc` | What gets dropped from MCParticle list |
| `lardataobj/Simulation/ParticleAncestryMap.h` | `sim::ParticleAncestryMap` for dropped particles |

---

## Verification

1. Build with `mrb build`
2. Run with existing simulation FCL
3. In Python: `import h5py; f = h5py.File('trackid_pid_map.h5'); list(f[frame_key].keys())`
   - Should now show: `track_ids`, `pids`, `mother_ids`, `processes`, `mother_pids`, `energy_fracs`
4. Sanity check:
   - Find a muon (PDG=13) → its daughters with process=Decay → one should be PDG=11 (Michel)
   - Find PDG=11 with process=eIoni and mother_pid=13/211 → Delta-ray
   - Confirm track count with MCParticle list is larger than SimChannel-only count
