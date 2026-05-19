#include "TrackIDPIDMap2h5.h"
#include "WireCellUtil/Configuration.h"
#include "WireCellUtil/NamedFactory.h"

#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "canvas/Utilities/InputTag.h"
#include "cetlib_except/exception.h"
#include "larsim/MCCheater/ParticleInventoryService.h"
#include "nusimdata/SimulationBase/MCParticle.h"

#include <algorithm>
#include <hdf5.h>
#include <map>
#include <string>

WIRECELL_FACTORY(TrackIDPIDMap2h5,
                 WireCell::AIML::TrackIDPIDMap2h5,
                 WireCell::INamed,
                 WireCell::IFrameFilter,
                 WireCell::IConfigurable)

using namespace WireCell;

// G4 process name → integer code mapping (following CellTree convention)
const std::map<std::string, int> AIML::TrackIDPIDMap2h5::s_process_map = {
  {"primary",         0},
  {"Decay",           1},
  {"eIoni",           2},
  {"muIoni",          3},
  {"eBrem",           4},
  {"compt",           5},
  {"phot",            6},
  {"conv",            7},
  {"hIoni",           8},
  {"nCapture",        9},
  {"muPairProd",     10},
  {"CoulombScat",    11},
  {"muBrems",        12},
  {"LowEnConversion",13},
  {"annihil",        14},
};

AIML::TrackIDPIDMap2h5::TrackIDPIDMap2h5()
  : Aux::Logger("TrackIDPIDMap2h5", "aiml")
  , m_simchannel_label("tpcrawdecoder:simpleSC")
  , m_particle_label("largeant")
  , m_output_file("trackid_pid_map.h5")
  , m_file(-1)
{}

AIML::TrackIDPIDMap2h5::~TrackIDPIDMap2h5()
{
  if (m_file >= 0) {
    H5Fclose(m_file);
    m_file = -1;
  }
}

Configuration AIML::TrackIDPIDMap2h5::default_configuration() const
{
  Configuration cfg;
  cfg["simchannel_label"] = m_simchannel_label;
  cfg["particle_label"] = m_particle_label;
  cfg["output_file"] = m_output_file;
  return cfg;
}

void AIML::TrackIDPIDMap2h5::configure(const Configuration& cfg)
{
  m_simchannel_label = get(cfg, "simchannel_label", m_simchannel_label);
  m_particle_label = get(cfg, "particle_label", m_particle_label);
  m_output_file = get(cfg, "output_file", m_output_file);

  clear_cache();
  if (m_file >= 0) {
    H5Fclose(m_file);
    m_file = -1;
  }
}

void AIML::TrackIDPIDMap2h5::visit(art::Event& event)
{
  log->info("TrackIDPIDMap2h5::visit called for event: {}", event.event());
  clear_cache();

  // --- Step 1: Read all MCParticles to get complete truth coverage ---
  // This covers ALL Geant4-tracked particles, not just TPC ionizers.
  if (!m_particle_label.empty()) {
    art::Handle<std::vector<simb::MCParticle>> particle_handle;
    if (event.getByLabel(art::InputTag{m_particle_label}, particle_handle)) {
      log->debug("TrackIDPIDMap2h5: Got {} MCParticles", particle_handle->size());
      for (auto const& particle : *particle_handle) {
        int tid = particle.TrackId();
        if (tid == 0) { continue; }
        m_trackid_to_pid[tid]      = particle.PdgCode();
        m_trackid_to_motherid[tid] = particle.Mother();
        // Map G4 process name to integer code; -1 for unknown
        auto it = s_process_map.find(particle.Process());
        m_trackid_to_process[tid] = (it != s_process_map.end()) ? it->second : -1;
        // Start/end positions from first and last trajectory points
        {
          size_t npts = particle.NumberTrajectoryPoints();
          const TLorentzVector& s4 = particle.Position(0);
          const TLorentzVector& e4 = particle.Position(npts - 1);
          m_trackid_to_start[tid] = {(float)s4.X(), (float)s4.Y(), (float)s4.Z(), (float)s4.T()};
          m_trackid_to_end[tid]   = {(float)e4.X(), (float)e4.Y(), (float)e4.Z(), (float)e4.T()};
        }
      }
      // Second pass: fill mother's PDG (requires full pid map built above)
      for (auto const& [tid, motherid] : m_trackid_to_motherid) {
        int mother_pdg = 0;
        if (motherid != 0) {
          auto mit = m_trackid_to_pid.find(motherid);
          if (mit != m_trackid_to_pid.end()) { mother_pdg = mit->second; }
        }
        m_trackid_to_motherpid[tid] = mother_pdg;
      }
    }
    else {
      log->warn("TrackIDPIDMap2h5 failed to fetch MCParticle with label '{}'", m_particle_label);
    }
  }

  // --- Step 2: Read SimChannel to accumulate energy fractions ---
  // Also ensures any trackID present only in SimChannel (e.g. shower secondaries
  // with modified IDs) is still recorded.
  if (!m_simchannel_label.empty()) {
    art::Handle<std::vector<sim::SimChannel>> handle;
    if (event.getByLabel(art::InputTag{m_simchannel_label}, handle)) {
      log->debug("TrackIDPIDMap2h5: Got {} SimChannels", handle->size());
      cache_simchannels(*handle);
      populate_trackid_pid_map();
    }
    else {
      log->warn("TrackIDPIDMap2h5 failed to fetch SimChannel with label '{}'", m_simchannel_label);
    }
  }

  log->info("TrackIDPIDMap2h5 cached {} track->pid entries", m_trackid_to_pid.size());
}

bool AIML::TrackIDPIDMap2h5::operator()(const input_pointer& in, output_pointer& out)
{
  out = in;
  if (!in) {
    log->debug("TrackIDPIDMap2h5: null input frame");
    return true;
  }

  log->info("TrackIDPIDMap2h5::operator() called for frame {}, have {} track->pid entries",
            in->ident(),
            m_trackid_to_pid.size());

  if (!m_trackid_to_pid.empty() || !m_simchnl_trackid_to_pid.empty()) {
    ensure_file();
    write_mapping(in->ident());
    write_mapping_simchnl(in->ident());
  }
  else {
    log->warn("TrackIDPIDMap2h5: no track ID to PID mapping to write for frame {} - was visit() called?",
              in->ident());
  }

  return true;
}

void AIML::TrackIDPIDMap2h5::cache_simchannels(const std::vector<sim::SimChannel>& simchs)
{
  m_simchannels = simchs;
}

void AIML::TrackIDPIDMap2h5::populate_trackid_pid_map()
{
  // Fill the SimChannel-based maps using ParticleInventoryService.
  // These cover only tracks that ionized in the TPC (including EM shower daughters
  // that are absent from the MCParticle list when keepEMShowerDaughters=false).
  if (m_simchannels.empty()) { return; }

  try {
    art::ServiceHandle<cheat::ParticleInventoryService> pi_serv;
    for (auto const& sc : m_simchannels) {
      for (auto const& tdc_entry : sc.TDCIDEMap()) {
        for (auto const& ide : tdc_entry.second) {
          const int track_id = ide.trackID;
          if (track_id == 0) { continue; }

          // Accumulate deposited energy [MeV] across all channels/TDCs
          m_simchnl_trackid_to_energy[track_id] += ide.energy;

          // Only look up once per trackID
          if (m_simchnl_trackid_to_pid.count(track_id)) { continue; }

          int pid = 0;
          int mother_id = 0;
          int process_code = -1;
          int mother_pdg = 0;
          try {
            auto const particle = pi_serv->TrackIdToParticle_P(track_id);
            if (particle) {
              pid = particle->PdgCode();
              mother_id = particle->Mother();
              auto it = s_process_map.find(particle->Process());
              process_code = (it != s_process_map.end()) ? it->second : -1;
              if (mother_id != 0) {
                auto const mother = pi_serv->TrackIdToParticle_P(mother_id);
                if (mother) { mother_pdg = mother->PdgCode(); }
              }
            }
          }
          catch (const cet::exception& ex) {
            log->debug("TrackIDPIDMap2h5: pi_serv lookup failed for track {}: {}",
                       track_id, ex.what());
          }
          m_simchnl_trackid_to_pid[track_id]       = pid;
          m_simchnl_trackid_to_motherid[track_id]  = mother_id;
          m_simchnl_trackid_to_process[track_id]   = process_code;
          m_simchnl_trackid_to_motherpid[track_id] = mother_pdg;
        }
      }
    }
    log->info("TrackIDPIDMap2h5: SimChannel-based map has {} entries", m_simchnl_trackid_to_pid.size());
  }
  catch (const cet::exception& ex) {
    log->warn("TrackIDPIDMap2h5: ParticleInventoryService unavailable: {}", ex.what());
  }
}

void AIML::TrackIDPIDMap2h5::clear_cache()
{
  m_simchannels.clear();
  // MCParticle-based maps
  m_trackid_to_pid.clear();
  m_trackid_to_motherid.clear();
  m_trackid_to_process.clear();
  m_trackid_to_motherpid.clear();
  m_trackid_to_start.clear();
  m_trackid_to_end.clear();
  // SimChannel-based maps
  m_simchnl_trackid_to_pid.clear();
  m_simchnl_trackid_to_motherid.clear();
  m_simchnl_trackid_to_process.clear();
  m_simchnl_trackid_to_motherpid.clear();
  m_simchnl_trackid_to_energy.clear();
}

void AIML::TrackIDPIDMap2h5::ensure_file()
{
  if (m_file >= 0) { return; }

  m_file = H5Fopen(m_output_file.c_str(), H5F_ACC_RDWR, H5P_DEFAULT);
  if (m_file < 0) {
    m_file = H5Fcreate(m_output_file.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
    if (m_file < 0) {
      log->error("TrackIDPIDMap2h5 failed to open output file {}", m_output_file);
    }
  }
}

// Helper: write a set of parallel int/float arrays into an HDF5 group
namespace {
  void h5_write_int(hid_t file, hid_t dataspace, hid_t lcpl,
                    const std::string& name, const std::vector<int>& data,
                    WireCell::Log::logptr_t log)
  {
    hid_t dset = H5Dcreate2(file, name.c_str(), H5T_NATIVE_INT, dataspace, lcpl, H5P_DEFAULT, H5P_DEFAULT);
    if (dset >= 0) {
      if (H5Dwrite(dset, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, data.data()) < 0)
        log->warn("TrackIDPIDMap2h5 failed to write {}", name);
      H5Dclose(dset);
    } else { log->warn("TrackIDPIDMap2h5 failed to create {}", name); }
  }
  void h5_write_float(hid_t file, hid_t dataspace, hid_t lcpl,
                      const std::string& name, const std::vector<float>& data,
                      WireCell::Log::logptr_t log)
  {
    hid_t dset = H5Dcreate2(file, name.c_str(), H5T_NATIVE_FLOAT, dataspace, lcpl, H5P_DEFAULT, H5P_DEFAULT);
    if (dset >= 0) {
      if (H5Dwrite(dset, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT, data.data()) < 0)
        log->warn("TrackIDPIDMap2h5 failed to write {}", name);
      H5Dclose(dset);
    } else { log->warn("TrackIDPIDMap2h5 failed to create {}", name); }
  }
  // Write a [N,4] float dataset (e.g. start/end XYZT positions)
  void h5_write_float_2d(hid_t file, hid_t lcpl,
                         const std::string& name,
                         const std::vector<std::array<float,4>>& data,
                         WireCell::Log::logptr_t log)
  {
    hsize_t dims[2] = {data.size(), 4};
    hid_t ds = H5Screate_simple(2, dims, nullptr);
    if (ds < 0) { log->warn("TrackIDPIDMap2h5 failed to create dataspace for {}", name); return; }
    hid_t dset = H5Dcreate2(file, name.c_str(), H5T_NATIVE_FLOAT, ds, lcpl, H5P_DEFAULT, H5P_DEFAULT);
    if (dset >= 0) {
      std::vector<float> flat;
      flat.reserve(data.size() * 4);
      for (auto const& a : data) { for (float v : a) flat.push_back(v); }
      if (H5Dwrite(dset, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT, flat.data()) < 0)
        log->warn("TrackIDPIDMap2h5 failed to write {}", name);
      H5Dclose(dset);
    } else { log->warn("TrackIDPIDMap2h5 failed to create {}", name); }
    H5Sclose(ds);
  }
} // anonymous namespace

void AIML::TrackIDPIDMap2h5::write_mapping(int frame_ident)
{
  // Writes MCParticle-based truth (all stored G4 tracks) under /<frame_ident>/mcpart/
  if (m_file < 0 || m_trackid_to_pid.empty()) { return; }

  std::vector<int> track_ids;
  track_ids.reserve(m_trackid_to_pid.size());
  for (auto const& entry : m_trackid_to_pid) { track_ids.push_back(entry.first); }
  std::sort(track_ids.begin(), track_ids.end());

  const size_t n = track_ids.size();
  std::vector<int> pids(n), mother_ids(n), processes(n), mother_pids(n);
  std::vector<std::array<float,4>> starts(n), ends(n);
  const std::array<float,4> zero4 = {0.f, 0.f, 0.f, 0.f};
  for (size_t i = 0; i < n; ++i) {
    int tid = track_ids[i];
    pids[i]        = m_trackid_to_pid.at(tid);
    mother_ids[i]  = m_trackid_to_motherid.count(tid)  ? m_trackid_to_motherid.at(tid)  : 0;
    processes[i]   = m_trackid_to_process.count(tid)   ? m_trackid_to_process.at(tid)   : -1;
    mother_pids[i] = m_trackid_to_motherpid.count(tid) ? m_trackid_to_motherpid.at(tid) : 0;
    starts[i]      = m_trackid_to_start.count(tid)     ? m_trackid_to_start.at(tid)     : zero4;
    ends[i]        = m_trackid_to_end.count(tid)       ? m_trackid_to_end.at(tid)       : zero4;
  }

  const hsize_t dims[1] = {n};
  const std::string grp = "/" + std::to_string(frame_ident) + "/mcpart";
  hid_t lcpl = H5Pcreate(H5P_LINK_CREATE);
  H5Pset_create_intermediate_group(lcpl, 1);
  hid_t dataspace = H5Screate_simple(1, dims, nullptr);
  if (dataspace < 0) { H5Pclose(lcpl); return; }

  h5_write_int(m_file, dataspace, lcpl, grp + "/track_ids",   track_ids,   log);
  h5_write_int(m_file, dataspace, lcpl, grp + "/pids",        pids,        log);
  h5_write_int(m_file, dataspace, lcpl, grp + "/mother_ids",  mother_ids,  log);
  h5_write_int(m_file, dataspace, lcpl, grp + "/processes",   processes,   log);
  h5_write_int(m_file, dataspace, lcpl, grp + "/mother_pids", mother_pids, log);

  H5Sclose(dataspace);
  // Write [N,4] start/end position datasets
  h5_write_float_2d(m_file, lcpl, grp + "/start_xyzts", starts, log);
  h5_write_float_2d(m_file, lcpl, grp + "/end_xyzts",   ends,   log);

  H5Pclose(lcpl);
  log->debug("TrackIDPIDMap2h5 wrote {} MCParticle entries to {}", n, grp);
}

void AIML::TrackIDPIDMap2h5::write_mapping_simchnl(int frame_ident)
{
  // Writes SimChannel-based truth (TPC ionizers, via ParticleInventoryService)
  // under /<frame_ident>/simchnl/
  if (m_file < 0 || m_simchnl_trackid_to_pid.empty()) { return; }

  std::vector<int> track_ids;
  track_ids.reserve(m_simchnl_trackid_to_pid.size());
  for (auto const& entry : m_simchnl_trackid_to_pid) { track_ids.push_back(entry.first); }
  std::sort(track_ids.begin(), track_ids.end());

  const size_t n = track_ids.size();
  std::vector<int>   pids(n), mother_ids(n), processes(n), mother_pids(n);
  std::vector<float> energies(n);
  for (size_t i = 0; i < n; ++i) {
    int tid = track_ids[i];
    pids[i]        = m_simchnl_trackid_to_pid.at(tid);
    mother_ids[i]  = m_simchnl_trackid_to_motherid.count(tid)  ? m_simchnl_trackid_to_motherid.at(tid)  : 0;
    processes[i]   = m_simchnl_trackid_to_process.count(tid)   ? m_simchnl_trackid_to_process.at(tid)   : -1;
    mother_pids[i] = m_simchnl_trackid_to_motherpid.count(tid) ? m_simchnl_trackid_to_motherpid.at(tid) : 0;
    energies[i]    = m_simchnl_trackid_to_energy.count(tid)    ? m_simchnl_trackid_to_energy.at(tid)    : 0.0f;
  }

  const hsize_t dims[1] = {n};
  const std::string grp = "/" + std::to_string(frame_ident) + "/simchnl";
  hid_t lcpl = H5Pcreate(H5P_LINK_CREATE);
  H5Pset_create_intermediate_group(lcpl, 1);
  hid_t dataspace = H5Screate_simple(1, dims, nullptr);
  if (dataspace < 0) { H5Pclose(lcpl); return; }

  h5_write_int  (m_file, dataspace, lcpl, grp + "/track_ids",   track_ids,   log);
  h5_write_int  (m_file, dataspace, lcpl, grp + "/pids",        pids,        log);
  h5_write_int  (m_file, dataspace, lcpl, grp + "/mother_ids",  mother_ids,  log);
  h5_write_int  (m_file, dataspace, lcpl, grp + "/processes",   processes,   log);
  h5_write_int  (m_file, dataspace, lcpl, grp + "/mother_pids", mother_pids, log);
  h5_write_float(m_file, dataspace, lcpl, grp + "/energies",    energies,    log);

  H5Sclose(dataspace);
  H5Pclose(lcpl);
  log->debug("TrackIDPIDMap2h5 wrote {} SimChannel entries to {}", n, grp);
}
