// snemo_gamma_tracking_efficiency_module.cc

// Ourselves:
#include <snemo_gamma_tracking_efficiency_module.h>

// Standard library:
#include <stdexcept>
#include <sstream>
#include <algorithm>

// Third party:
// - Boost:
#include <boost/fusion/iterator/next.hpp>

// - Bayeux/mctools
#include <mctools/utils.h>
#include <mctools/simulated_data.h>

// - Bayeux/datatools:
#include <datatools/clhep_units.h>
#include <datatools/service_manager.h>
// - Bayeux/mygsl
#include <mygsl/histogram_pool.h>
// - Bayeux/dpp
#include <dpp/histogram_service.h>
// - Bayeux/geomtools:
#include <geomtools/geometry_service.h>
#include <geomtools/manager.h>


// - Falaise
#include <snemo/processing/services.h>
#include <snemo/datamodels/data_model.h>
#include <snemo/datamodels/calibrated_data.h>
#include <snemo/datamodels/event_header.h>
#include <snemo/datamodels/particle_track.h>
#include <snemo/datamodels/particle_track_data.h>
#include <snemo/geometry/locator_plugin.h>
#include <snemo/geometry/calo_locator.h>
#include <snemo/geometry/xcalo_locator.h>
#include <snemo/geometry/gveto_locator.h>

namespace analysis {

  // Registration instantiation macro :
  DPP_MODULE_REGISTRATION_IMPLEMENT(snemo_gamma_tracking_efficiency_module,
                                    "analysis::snemo_gamma_tracking_efficiency");

  // Character separator between key for histogram dict.
  const char KEY_FIELD_SEPARATOR = '_';

  // Set the histogram pool used by the module :
  void snemo_gamma_tracking_efficiency_module::set_histogram_pool(mygsl::histogram_pool & pool_)
  {
    DT_THROW_IF(is_initialized(), std::logic_error,
                "Module '" << get_name() << "' is already initialized !");
    _histogram_pool_ = &pool_;
    return;
  }

  // Grab the histogram pool used by the module :
  mygsl::histogram_pool & snemo_gamma_tracking_efficiency_module::grab_histogram_pool()
  {
    DT_THROW_IF(! is_initialized(), std::logic_error,
                "Module '" << get_name() << "' is not initialized !");
    return *_histogram_pool_;
  }

  void snemo_gamma_tracking_efficiency_module::_set_defaults()
  {
    _key_fields_.clear ();

    _histogram_pool_ = 0;

    _efficiency_ = {0, 0, 0, 0, 0, 0, 0};

    _no_gt_efficiency_ = {0, 0};

    return;
  }

  // Initialization :
  void snemo_gamma_tracking_efficiency_module::initialize(const datatools::properties  & config_,
                                                          datatools::service_manager   & service_manager_,
                                                          dpp::module_handle_dict_type & module_dict_)
  {
    DT_THROW_IF(is_initialized(),
                std::logic_error,
                "Module '" << get_name() << "' is already initialized ! ");

    dpp::base_module::_common_initialize(config_);


    // Get the keys from 'Event Header' bank
    if (config_.has_key("key_fields"))
      {
        config_.fetch("key_fields", _key_fields_);
      }

    // Service label
    std::string histogram_label;
    if (config_.has_key("Histo_label"))
      {
        histogram_label = config_.fetch_string("Histo_label");
      }
    if (! _histogram_pool_)
      {
        DT_THROW_IF(histogram_label.empty(), std::logic_error,
                    "Module '" << get_name() << "' has no valid 'Histo_label' property !");

        DT_THROW_IF(! service_manager_.has(histogram_label) ||
                    ! service_manager_.is_a<dpp::histogram_service>(histogram_label),
                    std::logic_error,
                    "Module '" << get_name() << "' has no '" << histogram_label << "' service !");
        dpp::histogram_service & Histo
          = service_manager_.get<dpp::histogram_service>(histogram_label);
        set_histogram_pool(Histo.grab_pool());
        if (config_.has_key("Histo_output_files"))
          {
            std::vector<std::string> output_files;
            config_.fetch("Histo_output_files", output_files);
            for (size_t i = 0; i < output_files.size(); i++) {
              Histo.add_output_file(output_files[i]);
            }
          }
        if (config_.has_key("Histo_input_file"))
          {
            const std::string input_file = config_.fetch_string("Histo_input_file");
            Histo.load_from_boost_file(input_file);
          }
        if (config_.has_key("Histo_template_files"))
          {
            std::vector<std::string> template_files;
            config_.fetch("Histo_template_files", template_files);
            for (size_t i = 0; i < template_files.size(); i++) {
              Histo.grab_pool().load(template_files[i]);
            }
          }

        // Geometry manager :
        std::string geo_label = snemo::processing::service_info::default_geometry_service_label();
        if (config_.has_key("Geo_label")) {
          geo_label = config_.fetch_string("Geo_label");
        }
        DT_THROW_IF (geo_label.empty(), std::logic_error,
                     "Module '" << get_name() << "' has no valid '" << "Geo_label" << "' property !");
        DT_THROW_IF (! service_manager_.has(geo_label) ||
                     ! service_manager_.is_a<geomtools::geometry_service>(geo_label),
                     std::logic_error,
                     "Module '" << get_name() << "' has no '" << geo_label << "' service !");
        geomtools::geometry_service & Geo
          = service_manager_.get<geomtools::geometry_service>(geo_label);

        // Get geometry locator plugin
        const geomtools::manager & geo_mgr = Geo.get_geom_manager();
        std::string locator_plugin_name;
        if (config_.has_key ("locator_plugin_name"))
          {
            locator_plugin_name = config_.fetch_string ("locator_plugin_name");
          }
        else
          {
            // If no locator plugin name is set, then search for the first one
            const geomtools::manager::plugins_dict_type & plugins = geo_mgr.get_plugins ();
            for (geomtools::manager::plugins_dict_type::const_iterator ip = plugins.begin ();
                 ip != plugins.end ();
                 ip++) {
              const std::string & plugin_name = ip->first;
              if (geo_mgr.is_plugin_a<snemo::geometry::locator_plugin> (plugin_name)) {
                DT_LOG_DEBUG (get_logging_priority (), "Find locator plugin with name = " << plugin_name);
                locator_plugin_name = plugin_name;
                break;
              }
            }
          }
        // Access to a given plugin by name and type :
        DT_THROW_IF (! geo_mgr.has_plugin (locator_plugin_name) ||
                     ! geo_mgr.is_plugin_a<snemo::geometry::locator_plugin> (locator_plugin_name),
                     std::logic_error,
                     "Found no locator plugin named '" << locator_plugin_name << "'");
        _locator_plugin_ = &geo_mgr.get_plugin<snemo::geometry::locator_plugin> (locator_plugin_name);

        // Tag the module as initialized :
        _set_initialized(true);
        return;
      }
  }
  // Reset :
  void snemo_gamma_tracking_efficiency_module::reset()
  {
    DT_THROW_IF(! is_initialized(), std::logic_error,
                "Module '" << get_name() << "' is not initialized !");

    // Present results
    DT_LOG_NOTICE(get_logging_priority(),
                  "Number of gammas well reconstructed = " << _efficiency_.ngood << " / " << _efficiency_.ntotal
                  << " ( " << _efficiency_.ngood/(double)_efficiency_.ntotal*100 << " %)");
    DT_LOG_NOTICE(get_logging_priority(),
                  "Number of gammas missed = " << _efficiency_.nmiss << " / " << _efficiency_.nevent
                  << " ( " << _efficiency_.nmiss/(double)_efficiency_.nevent*100 << " %)");

    DT_LOG_NOTICE(get_logging_priority(),
                  "Number of events successfully reconstructed = " << _efficiency_.ngood_event << " / " << _efficiency_.nevent
                  << " ( " << _efficiency_.ngood_event/(double)_efficiency_.nevent*100 << " %)");
    DT_LOG_NOTICE(get_logging_priority(),
                  "Number of events with gammas successfully reconstructed = " << _efficiency_.ngood_event << " / " << _efficiency_.nevent_gammas
                  << " ( " << _efficiency_.ngood_event/(double)_efficiency_.nevent_gammas*100 << " %)");

    DT_LOG_WARNING(get_logging_priority(),
                   "Number of events with gammas successfully reconstructed = " << _efficiency_.ngood_event << " / " << _efficiency_.nevent_gammas
                   << " ( " << _efficiency_.ngood_event/(double)_efficiency_.nevent_gammas*100 << " %)");

    DT_LOG_WARNING(get_logging_priority(),
                   "Number of events with gammas successfully clustered = " << _no_gt_efficiency_.no_gt_ngood_event << " / " << _no_gt_efficiency_.no_gt_nevent_gammas
                   << " ( " << _no_gt_efficiency_.no_gt_ngood_event/(double)_no_gt_efficiency_.no_gt_nevent_gammas*100 << " %)");

    // Tag the module as un-initialized :
    _set_initialized(false);
    _set_defaults();
    return;
  }

  // Constructor :
  snemo_gamma_tracking_efficiency_module::snemo_gamma_tracking_efficiency_module(datatools::logger::priority logging_priority_)
    : dpp::base_module(logging_priority_)
  {
    _set_defaults();
    return;
  }

  // Destructor :
  snemo_gamma_tracking_efficiency_module::~snemo_gamma_tracking_efficiency_module()
  {
    if (is_initialized()) snemo_gamma_tracking_efficiency_module::reset();
    return;
  }

  // Explore the cluster
  void snemo_gamma_tracking_efficiency_module::get_new_neighbours(geomtools::geom_id gid,
                                                                  const snemo::datamodel::calibrated_data::calorimeter_hit_collection_type & cch,
                                                                  std::vector<geomtools::geom_id>  & ccl,
                                                                  std::vector<geomtools::geom_id>  & a_cluster)
  {
    if(std::find(ccl.begin(), ccl.end(),gid)==ccl.end())
      ccl.push_back(gid);
    else
      return;

    std::vector<geomtools::geom_id>  the_neighbours = {};
    std::vector<geomtools::geom_id>  the_calib_neighbours = {};

    const snemo::geometry::calo_locator & calo_locator
      = _locator_plugin_->get_calo_locator();
    const snemo::geometry::xcalo_locator & xcalo_locator
      = _locator_plugin_->get_xcalo_locator();
    const snemo::geometry::gveto_locator & gveto_locator
      = _locator_plugin_->get_gveto_locator();

    if (calo_locator.is_calo_block_in_current_module(gid))
      calo_locator.get_neighbours_ids(gid, the_neighbours, snemo::geometry::utils::NEIGHBOUR_FIRST);

    if (xcalo_locator.is_calo_block_in_current_module(gid))
      xcalo_locator.get_neighbours_ids(gid, the_neighbours, snemo::geometry::utils::NEIGHBOUR_FIRST);

    if (gveto_locator.is_calo_block_in_current_module(gid))
      gveto_locator.get_neighbours_ids(gid, the_neighbours, snemo::geometry::utils::NEIGHBOUR_FIRST);

    for (auto ineighbour : the_neighbours)
      if (std::find_if(cch.begin(), cch.end(), [ineighbour] (auto icalo)
                       {return ineighbour == icalo.get().get_geom_id();}) != cch.end())
        if(std::find(ccl.begin(), ccl.end(), ineighbour)==ccl.end()) {
          the_calib_neighbours.push_back(ineighbour);
          ccl.push_back(ineighbour);
          a_cluster.push_back(ineighbour);
        }

    for(auto i_calib_neighbour : the_calib_neighbours)
      get_new_neighbours(i_calib_neighbour, cch, ccl, a_cluster);
  }

  // Pre processing for cluster identification
  void snemo_gamma_tracking_efficiency_module::_pre_process_clustering(const datatools::things & data_record_,
                                                                       gamma_dict_type & clustered_gammas_)
  {
    // // Check if some 'calibrated_data' are available in the data model:
    // const std::string cd_label = snemo::datamodel::data_info::default_calibrated_data_label();
    // if (! data_record_.has(cd_label)) {
    //   DT_LOG_ERROR(get_logging_priority(), "Missing calibrated data to be processed !");
    // }

    // // Get the 'calibrated_data' entry from the data model :
    // const snemo::datamodel::calibrated_data & cd
    //   = data_record_.get<snemo::datamodel::calibrated_data>(cd_label);

    // Check if some 'particle_track_data' are available in the data model:
    const std::string ptd_label = snemo::datamodel::data_info::default_particle_track_data_label();
    if (! data_record_.has(ptd_label)) {
      DT_LOG_ERROR(get_logging_priority(), "Missing particle track data to be processed !");
    }

    // Get the 'particle_track_data' entry from the data model :
    const snemo::datamodel::particle_track_data & ptd
      = data_record_.get<snemo::datamodel::particle_track_data>(ptd_label);

    snemo::datamodel::particle_track_data::particle_collection_type gamma_particles;
    ptd.fetch_particles(gamma_particles, snemo::datamodel::particle_track::NEUTRAL);

    // retrieve only hits from gammas
    snemo::datamodel::calibrated_calorimeter_hit::collection_type cch;

    for(auto igamma : gamma_particles) {
      snemo::datamodel::particle_track a_gamma = igamma.grab();
      cch.insert(cch.end(),a_gamma.get_associated_calorimeter_hits().begin(),a_gamma.get_associated_calorimeter_hits().end());
    }

    // std::cout << " cch size " << cch.size() << std::endl;

    std::vector<geomtools::geom_id>  ccl = {};

    size_t number_of_clusters = 0;

    std::vector<std::vector<geomtools::geom_id> >  the_reconstructed_clusters;

    for (auto icalo : cch) {

      const geomtools::geom_id & gid = icalo.get().get_geom_id();

      std::vector<geomtools::geom_id> a_cluster = {};
      a_cluster.push_back(gid);

      if(std::find(ccl.begin(), ccl.end(),gid)!=ccl.end())
        continue;

      get_new_neighbours(gid, cch, ccl, a_cluster);

      the_reconstructed_clusters.push_back(a_cluster);

      number_of_clusters++;
    }

    std::vector<std::map<double, geomtools::geom_id> >  the_ordered_reconstructed_clusters;

    for(auto icluster : the_reconstructed_clusters)
      {
        std::map<double,geomtools::geom_id> a_cluster = {};

       for(auto igid : icluster)
          for(auto ihit : cch)
            if(igid == ihit.get().get_geom_id())
                a_cluster.insert( std::pair<double,geomtools::geom_id >(ihit.get().get_time(),igid) );

       the_ordered_reconstructed_clusters.push_back(a_cluster);
      }

    std::sort(the_ordered_reconstructed_clusters.begin(), the_ordered_reconstructed_clusters.end());

    int track_id = 0;

    for(auto icluster : the_ordered_reconstructed_clusters)
      {
        track_id++;

        if(icluster.size() < 2)
          {
            for (auto ipair : icluster)
              clustered_gammas_[track_id].insert(ipair.second);
            continue;
          }

        double t0 = 0; // not ideal
        double t1 = 0;

        for (auto ipair : icluster)
          {
            t0 = t1;
            t1 = ipair.first;
            // std::cout << " " << ipair.first  << "   " << std::next(&ipair)->first << std::endl;
            // std::cout << " " << t0  << "   " << t1 << std::endl;

            if(t0!=0 && t1!=0 && t1-t0 > 2.5 /*ns*/)
              {
                number_of_clusters++;
                track_id++;
              }

            clustered_gammas_[track_id].insert(ipair.second);
          }
      }

    // Check if some 'event_header' are available in the data model:
    const std::string eh_label = snemo::datamodel::data_info::default_event_header_label();
    if (! data_record_.has(eh_label)) {
      DT_LOG_ERROR(get_logging_priority(), "Missing event header info !");
    }

    // Get the 'event_header' entry from the data model :
    const snemo::datamodel::event_header & eh
      = data_record_.get<snemo::datamodel::event_header>(eh_label);

    // if(number_of_clusters == 4)
    //   std::cout << eh.get_id() << std::endl;

    // Build unique key for histogram map:
    std::ostringstream key;
    key << "number_of_gamma_calos";

    // Getting histogram pool
    mygsl::histogram_pool & a_pool = grab_histogram_pool();

    if (! a_pool.has(key.str()))
      {
        mygsl::histogram_1d & h = a_pool.add_1d(key.str(), "", "number_of_gamma_calos");
        datatools::properties hconfig;
        hconfig.store_string("mode", "mimic");
        hconfig.store_string("mimic.histogram_1d", "number_of_calos_template");
        mygsl::histogram_pool::init_histo_1d(h, hconfig, &a_pool);
      }

    // Getting the current histogram
    mygsl::histogram_1d & a_histo_gamma_calos = a_pool.grab_1d(key.str ());
    a_histo_gamma_calos.fill(cch.size());

    key.str("");
    key.clear();

    // // Build unique key for histogram map:
    // std::ostringstream key;
    key << "number_of_gamma_clusters";

    // // Getting histogram pool
    // mygsl::histogram_pool & a_pool = grab_histogram_pool();

    if (! a_pool.has(key.str()))
      {
        mygsl::histogram_1d & h = a_pool.add_1d(key.str(), "", "number_of_gamma_clusters");
        datatools::properties hconfig;
        hconfig.store_string("mode", "mimic");
        hconfig.store_string("mimic.histogram_1d", "number_of_calos_template");
        mygsl::histogram_pool::init_histo_1d(h, hconfig, &a_pool);
      }

    // Getting the current histogram
    mygsl::histogram_1d & a_histo_gamma_clusters = a_pool.grab_1d(key.str ());
    a_histo_gamma_clusters.fill(number_of_clusters);

        key.str("");
    key.clear();

    // // Build unique key for histogram map:
    // std::ostringstream key;
    key << "clusters_size";

    // // Getting histogram pool
    // mygsl::histogram_pool & a_pool = grab_histogram_pool();

    if (! a_pool.has(key.str()))
      {
        mygsl::histogram_1d & h = a_pool.add_1d(key.str(), "", "clusters_size");
        datatools::properties hconfig;
        hconfig.store_string("mode", "mimic");
        hconfig.store_string("mimic.histogram_1d", "number_of_calos_template");
        mygsl::histogram_pool::init_histo_1d(h, hconfig, &a_pool);
      }

    // Getting the current histogram
    mygsl::histogram_1d & a_histo_clusters_size = a_pool.grab_1d(key.str ());

    for (auto igamma : clustered_gammas_)
        a_histo_clusters_size.fill(igamma.second.size());

    // _no_gt_efficiency_.no_gt_ngood_event++;

  }

// Processing :
dpp::base_module::process_status snemo_gamma_tracking_efficiency_module::process(datatools::things & data_record_)
{
  DT_LOG_TRACE(get_logging_priority(), "Entering...");
  DT_THROW_IF(! is_initialized(), std::logic_error,
              "Module '" << get_name() << "' is not initialized !");

   // std::cout << " ---------------------------------------------------------------------------------- " << std::endl;

  gamma_dict_type clustered_gammas;
  {
    _pre_process_clustering(data_record_, clustered_gammas);
  }

  gamma_dict_type simulated_gammas;
  {
    const process_status status = _process_simulated_gammas(data_record_, simulated_gammas);
    if (status != dpp::base_module::PROCESS_OK) {
      DT_LOG_ERROR(get_logging_priority(), "Processing of simulated data fails !");
      return status;
    }
  }

  gamma_dict_type reconstructed_gammas;
  {
    const process_status status = _process_reconstructed_gammas(data_record_, reconstructed_gammas);
    if (status != dpp::base_module::PROCESS_OK) {
      DT_LOG_ERROR(get_logging_priority(), "Processing of particle track data fails !");
      return status;
    }
    // return dpp::base_module::PROCESS_OK;
  }

  _compare_sequences(simulated_gammas, reconstructed_gammas);

  _compare_sequences_cluster(simulated_gammas, clustered_gammas);

  // const process_status status = _compute_gamma_track_length(data_record_);
  // if (status != dpp::base_module::PROCESS_OK) {
  //   DT_LOG_ERROR(get_logging_priority(), "Computing the gamma track length fails !");
  //   return status;
  // }

  DT_LOG_TRACE(get_logging_priority(), "Exiting.");
    return dpp::base_module::PROCESS_SUCCESS;

}

dpp::base_module::process_status snemo_gamma_tracking_efficiency_module::_process_simulated_gammas(const datatools::things & data_record_,
                                                                                                   gamma_dict_type & simulated_gammas_)
{
  // Check if some 'simulated_data' are available in the data model:
  const std::string sd_label = snemo::datamodel::data_info::default_simulated_data_label();
  if (! data_record_.has(sd_label)) {
    DT_LOG_ERROR(get_logging_priority(), "Missing simulated data to be processed !");
    return dpp::base_module::PROCESS_ERROR;
  }
  // Get the 'simulated_data' entry from the data model :
  const mctools::simulated_data & sd = data_record_.get<mctools::simulated_data>(sd_label);

  DT_LOG_DEBUG(get_logging_priority(), "Simulated data : ");
  if (get_logging_priority() >= datatools::logger::PRIO_DEBUG) sd.tree_dump();

  // Get total number of gammas simulated
  _efficiency_.ngamma = 0;
  for (auto i : sd.get_primary_event().get_particles()) {
    if (i.is_gamma()) _efficiency_.ngamma++;
  }

  // Check if some 'calibrated_data' are available in the data model:
  const std::string cd_label = snemo::datamodel::data_info::default_calibrated_data_label();
  if (! data_record_.has(cd_label)) {
    DT_LOG_ERROR(get_logging_priority(), "Missing calibrated data to be processed !");
    return dpp::base_module::PROCESS_ERROR;
  }

  // Get the 'calibrated_data' entry from the data model :
  const snemo::datamodel::calibrated_data & cd
    = data_record_.get<snemo::datamodel::calibrated_data>(cd_label);

  DT_LOG_DEBUG(get_logging_priority(), "Calibrated data : ");
  if (get_logging_priority() >= datatools::logger::PRIO_DEBUG) cd.tree_dump();

  // Stop proccess if no calibrated calorimeters
  if (! cd.has_calibrated_calorimeter_hits())
    return dpp::base_module::PROCESS_STOP;

  const snemo::datamodel::calibrated_data::calorimeter_hit_collection_type & cch
    = cd.calibrated_calorimeter_hits();

  // Fetch simulated step hits from calorimeter blocks
  const std::string hit_label = "__visu.tracks.calo";
  if (! sd.has_step_hits(hit_label)) return dpp::base_module::PROCESS_STOP;
  const mctools::simulated_data::hit_handle_collection_type & hit_collection
    = sd.get_step_hits(hit_label);
  if (hit_collection.empty()) {
    DT_LOG_DEBUG(get_logging_priority(), "No simulated calorimeter hits");
    return dpp::base_module::PROCESS_STOP;
  }

  std::set<geomtools::geom_id> already_gids;

  for (auto ihit : hit_collection) {
    const mctools::base_step_hit & a_hit = ihit.get();
    const datatools::properties & a_aux = a_hit.get_auxiliaries();

    int track_id = -1;
    if (a_aux.has_key(mctools::track_utils::TRACK_ID_KEY)) {
      track_id = a_aux.fetch_integer(mctools::track_utils::TRACK_ID_KEY);
    }
    if (a_aux.has_key(mctools::track_utils::PARENT_TRACK_ID_KEY)) {
      track_id = a_aux.fetch_integer(mctools::track_utils::PARENT_TRACK_ID_KEY);
    }
    DT_THROW_IF(track_id == -1, std::logic_error, "Missing primary track id !");
    if (track_id == 0) continue; // From a primary particles

    // Check if calorimeter has been calibrated
    const geomtools::geom_id & gid = a_hit.get_geom_id();
    if (
        std::find_if(cch.begin(), cch.end(), [gid] (auto hit_)
                     {
                       return gid == hit_.get().get_geom_id();
                     }
                     ) == cch.end()
        ) continue;


    // if(simulated_gammas_[track_id].size() != 0 )
    //   {
    //     // DT_LOG_WARNING(get_logging_priority(), std::endl << "track id already added !" << std::endl);
    //     continue;
    //   }

    // Gid already attributed to a gamma
    if (already_gids.count(gid) != 0)
      {
        // DT_LOG_WARNING(get_logging_priority(), std::endl << "GID already found !" << already_gids.count(gid) << std::endl);
        continue;
      }

    if (track_id > (int)_efficiency_.ngamma + 1) //continue; // Not from a primary particles // Hack : removes around 10% of the stat
      {
        DT_LOG_WARNING(get_logging_priority(), "Secondary particle triggering new calo "<< _efficiency_.ngamma);
        return dpp::base_module::PROCESS_STOP;
      }

    already_gids.insert(gid);

    // DT_LOG_WARNING(get_logging_priority(), std::endl << "Inserting gid " << already_gids.count(gid));

    simulated_gammas_[track_id].insert(gid);

    // Build unique key for histogram map:
    std::ostringstream key;
    key << "total_number_of_calos";

    // Getting histogram pool
    mygsl::histogram_pool & a_pool = grab_histogram_pool();

    if (! a_pool.has(key.str()))
      {
        mygsl::histogram_1d & h = a_pool.add_1d(key.str(), "", "number_of_calos");
        datatools::properties hconfig;
        hconfig.store_string("mode", "mimic");
        hconfig.store_string("mimic.histogram_1d", "number_of_calos_template");
        mygsl::histogram_pool::init_histo_1d(h, hconfig, &a_pool);
      }

    // Getting the current histogram
    mygsl::histogram_1d & a_histo = a_pool.grab_1d(key.str ());
    a_histo.fill(cch.size());

    // std::cout << " simulated cch size " << cch.size() <<std::endl;


  }

  return dpp::base_module::PROCESS_OK;
}

dpp::base_module::process_status snemo_gamma_tracking_efficiency_module::_compute_gamma_track_length(const datatools::things & data_record_)
{

  // Check if some 'simulated_data' are available in the data model:
  const std::string sd_label = snemo::datamodel::data_info::default_simulated_data_label();
  if (! data_record_.has(sd_label)) {
    DT_LOG_ERROR(get_logging_priority(), "Missing simulated data to be processed !");
    return dpp::base_module::PROCESS_ERROR;
  }
  // Get the 'simulated_data' entry from the data model :
  const mctools::simulated_data & sd = data_record_.get<mctools::simulated_data>(sd_label);

  double simu_gamma_track_length = 0.;
  for (auto i : sd.get_primary_event().get_particles()) {

    if (!i.is_gamma()) continue;

    for(unsigned int i = 0; i<1; i++) //if more than one label
      {
        // Fetch simulated step hits from calorimeter blocks
        // const std::string hit_label = "__visu.tracks";
        const std::string hit_label = "__visu.tracks.calo";
        if (! sd.has_step_hits(hit_label)) return dpp::base_module::PROCESS_STOP;
        const mctools::simulated_data::hit_handle_collection_type & hit_collection
          = sd.get_step_hits(hit_label);
        if (hit_collection.empty())
          {
            DT_LOG_DEBUG(get_logging_priority(), "No simulated calorimeter hits");
            return dpp::base_module::PROCESS_STOP;
          }

        for (auto ihit : hit_collection)
          {
            const mctools::base_step_hit & a_hit = ihit.get();
            const datatools::properties & a_aux = a_hit.get_auxiliaries();

            const geomtools::vector_3d & a_hit_start = a_hit.get_position_start();
            const geomtools::vector_3d & a_hit_stop = a_hit.get_position_stop();

            simu_gamma_track_length += (a_hit_stop-a_hit_start).mag();
          }
      }
    std::cout<< "------- Simulated gamma track length : " << simu_gamma_track_length << std::endl;
  }

  // Check if some 'particle_track_data' are available in the data model:
  const std::string ptd_label = snemo::datamodel::data_info::default_particle_track_data_label();
  if (! data_record_.has(ptd_label)) {
    DT_LOG_ERROR(get_logging_priority(), "Missing particle track data to be processed !");
    return dpp::base_module::PROCESS_ERROR;
  }

  // Get the 'particle_track_data' entry from the data model :
  const snemo::datamodel::particle_track_data & ptd
    = data_record_.get<snemo::datamodel::particle_track_data>(ptd_label);

  DT_LOG_DEBUG(get_logging_priority(), "Particle track data : ");
  if (get_logging_priority() >= datatools::logger::PRIO_DEBUG) ptd.tree_dump();

  snemo::datamodel::particle_track_data::particle_collection_type gammas;
  const size_t ngammas = ptd.fetch_particles(gammas, snemo::datamodel::particle_track::NEUTRAL);

  if (ngammas == 0) return dpp::base_module::PROCESS_STOP;

  double reco_gamma_track_length = 0;

  for (auto igamma : gammas) {


    const snemo::datamodel::particle_track::vertex_collection_type the_vertices = igamma.get().get_vertices();

    unsigned int count_vtx = 0;

    for (auto ivtx : the_vertices)
      {
        // count_vtx++;

        // if (count_vtx == the_vertices.size()-1)
        //   break;

        const geomtools::blur_spot & a_spot = ivtx.get();
        const geomtools::blur_spot & next_spot = std::next(&ivtx)->get();

        reco_gamma_track_length += (next_spot.get_position() - a_spot.get_position()).mag();

      }

    std::cout<< "------- Reconstructed gamma track length : " << reco_gamma_track_length << std::endl;


    // Check if some 'calibrated_data' are available in the data model:
    const std::string cd_label = snemo::datamodel::data_info::default_calibrated_data_label();
    if (! data_record_.has(cd_label)) {
      DT_LOG_ERROR(get_logging_priority(), "Missing calibrated data to be processed !");
      return dpp::base_module::PROCESS_ERROR;
    }
    // Get the 'calibrated_data' entry from the data model :
    const snemo::datamodel::calibrated_data & cd
      = data_record_.get<snemo::datamodel::calibrated_data>(cd_label);

    // Stop proccess if no calibrated calorimeters
    if (! cd.has_calibrated_calorimeter_hits()) return dpp::base_module::PROCESS_STOP;

    const snemo::datamodel::calibrated_data::calorimeter_hit_collection_type & cch
      = cd.calibrated_calorimeter_hits();

    // Build unique key for histogram map:
    std::ostringstream key;
    key << cch.size();
    key <<"calos_";
    key << "delta_L";

    if(reco_gamma_track_length < 400)
      key << "_cluster";
    else
      key << "not_cluster";

    if(cch.size()==2 && (simu_gamma_track_length - reco_gamma_track_length ) < -800.)
      std::cout << std::endl << "Check it out " << std::endl << std::endl;

    // Getting histogram pool
    mygsl::histogram_pool & a_pool = grab_histogram_pool();

    if (! a_pool.has(key.str()))
      {
        mygsl::histogram_1d & h = a_pool.add_1d(key.str(), "", "delta_L");
        datatools::properties hconfig;
        hconfig.store_string("mode", "mimic");
        hconfig.store_string("mimic.histogram_1d", "delta_L_template");
        mygsl::histogram_pool::init_histo_1d(h, hconfig, &a_pool);
      }

    // Getting the current histogram
    mygsl::histogram_1d & a_histo = a_pool.grab_1d(key.str ());
    a_histo.fill((simu_gamma_track_length - reco_gamma_track_length )*CLHEP::keV);

  }

  return dpp::base_module::PROCESS_OK;
}

dpp::base_module::process_status snemo_gamma_tracking_efficiency_module::_process_reconstructed_gammas(const datatools::things & data_record_,
                                                                                                       gamma_dict_type & reconstructed_gammas_)
{
  // Check if some 'particle_track_data' are available in the data model:
  const std::string ptd_label = snemo::datamodel::data_info::default_particle_track_data_label();
  if (! data_record_.has(ptd_label)) {
    DT_LOG_ERROR(get_logging_priority(), "Missing particle track data to be processed !");
    return dpp::base_module::PROCESS_ERROR;
  }

  // Get the 'particle_track_data' entry from the data model :
  const snemo::datamodel::particle_track_data & ptd
    = data_record_.get<snemo::datamodel::particle_track_data>(ptd_label);

  DT_LOG_DEBUG(get_logging_priority(), "Particle track data : ");
  if (get_logging_priority() >= datatools::logger::PRIO_DEBUG) ptd.tree_dump();

  snemo::datamodel::particle_track_data::particle_collection_type gammas;
  const size_t ngammas = ptd.fetch_particles(gammas, snemo::datamodel::particle_track::NEUTRAL);
  if (ngammas == 0) return dpp::base_module::PROCESS_STOP;

  DT_LOG_DEBUG(get_logging_priority(), std::endl << "Number of gammas : " << ngammas << std::endl);

  for (auto igamma : gammas) {
    for (auto icalo : igamma.get().get_associated_calorimeter_hits()) {
      const geomtools::geom_id & gid = icalo.get().get_geom_id();
      reconstructed_gammas_[igamma.get().get_track_id()].insert(gid);
    }
  }

  // Build unique key for histogram map:
  std::ostringstream key;
  key << "number_of_gammas";

  // Getting histogram pool
  mygsl::histogram_pool & a_pool = grab_histogram_pool();

  if (! a_pool.has(key.str()))
    {
      mygsl::histogram_1d & h = a_pool.add_1d(key.str(), "", "number_of_calos");
      datatools::properties hconfig;
      hconfig.store_string("mode", "mimic");
      hconfig.store_string("mimic.histogram_1d", "number_of_calos_template");
      mygsl::histogram_pool::init_histo_1d(h, hconfig, &a_pool);
    }

  // Getting the current histogram
  mygsl::histogram_1d & a_histo = a_pool.grab_1d(key.str ());
  a_histo.fill(ngammas);

  key.str("");
  key.clear();

  key << "total_gamma_energy";

  double total_gamma_energy = 0;
  double E_1 = 0;
  double E_2 = 0;
  double E_3 = 0;
  double E_min = 0;
  double E_mid = 0;
  double E_max = 0;
  double ncalos_1 = 0;
  double ncalos_2 = 0;
  double ncalos_3 = 0;
  double ncalos_min = 0;
  double ncalos_mid = 0;
  double ncalos_max = 0;

  int count_gamma = 0;
  for (auto igamma : gammas) {
    count_gamma++;
    // std::cout << "igamma : " << igamma << std::endl;
    for (auto icalo : igamma.get().get_associated_calorimeter_hits()) {
      total_gamma_energy += icalo.get().get_energy();
      if(count_gamma==1)
        {
          E_1 += icalo.get().get_energy();
          ncalos_1 = igamma.get().get_associated_calorimeter_hits().size();
        }
      if(count_gamma==2)
        {
          E_2 += icalo.get().get_energy();
          ncalos_2 = igamma.get().get_associated_calorimeter_hits().size();
        }
      if(count_gamma==3)
        {
          E_3 += icalo.get().get_energy();
          ncalos_3 = igamma.get().get_associated_calorimeter_hits().size();
        }
    }
  }

  // std::cout << "E1 : " << E_1 << std::endl;
  // std::cout << "E2 : " << E_2 << std::endl;
  // std::cout << "E3 : " << E_3 << std::endl;

  if(ngammas == 3)
    {
      E_min = std::min(std::min(E_1,E_2),E_3);
      E_mid = std::min(std::max(E_1,E_2),E_3);
      E_max = std::max(std::max(E_1,E_2),E_3);
    }

  if (! a_pool.has(key.str()))
    {
      mygsl::histogram_1d & h = a_pool.add_1d(key.str(), "", "total_gamma_energy");
      datatools::properties hconfig;
      hconfig.store_string("mode", "mimic");
      hconfig.store_string("mimic.histogram_1d", "energy_template");
      mygsl::histogram_pool::init_histo_1d(h, hconfig, &a_pool);
    }

  // Getting the current histogram
  mygsl::histogram_1d & a_histo_total_gamma_energy = a_pool.grab_1d(key.str ());
  a_histo_total_gamma_energy.fill(total_gamma_energy);

  if(E_min != 0)
    {
      key.str("");
      key.clear();

      if(E_min == E_1)
        key << ncalos_1;
      if(E_min == E_2)
        key << ncalos_2;
      if(E_min == E_3)
        key << ncalos_3;

      key << "_gamma_energy_min";

      if (! a_pool.has(key.str()))
        {
          mygsl::histogram_1d & h = a_pool.add_1d(key.str(), "", "gamma_energy_min");
          datatools::properties hconfig;
          hconfig.store_string("mode", "mimic");
          hconfig.store_string("mimic.histogram_1d", "energy_template");
          mygsl::histogram_pool::init_histo_1d(h, hconfig, &a_pool);
        }

      // Getting the current histogram
      mygsl::histogram_1d & a_histo_gamma_energy_min = a_pool.grab_1d(key.str ());
      a_histo_gamma_energy_min.fill(E_min);
    }

  if(E_mid != 0)
    {
      key.str("");
      key.clear();

      if(E_min == E_1)
        key << ncalos_1;
      if(E_min == E_2)
        key << ncalos_2;
      if(E_min == E_3)
        key << ncalos_3;

      key << "_gamma_energy_mid";

      if (! a_pool.has(key.str()))
        {
          mygsl::histogram_1d & h = a_pool.add_1d(key.str(), "", "gamma_energy_mid");
          datatools::properties hconfig;
          hconfig.store_string("mode", "mimic");
          hconfig.store_string("mimic.histogram_1d", "energy_template");
          mygsl::histogram_pool::init_histo_1d(h, hconfig, &a_pool);
        }

      // Getting the current histogram
      mygsl::histogram_1d & a_histo_gamma_energy_mid = a_pool.grab_1d(key.str ());
      a_histo_gamma_energy_mid.fill(E_mid);
    }

  if(E_max != 0)
    {
      key.str("");
      key.clear();

      if(E_min == E_1)
        key << ncalos_1;
      if(E_min == E_2)
        key << ncalos_2;
      if(E_min == E_3)
        key << ncalos_3;

      key << "_gamma_energy_max";

      if (! a_pool.has(key.str()))
        {
          mygsl::histogram_1d & h = a_pool.add_1d(key.str(), "", "gamma_energy_max");
          datatools::properties hconfig;
          hconfig.store_string("mode", "mimic");
          hconfig.store_string("mimic.histogram_1d", "energy_template");
          mygsl::histogram_pool::init_histo_1d(h, hconfig, &a_pool);
        }

      // Getting the current histogram
      mygsl::histogram_1d & a_histo_gamma_energy_max = a_pool.grab_1d(key.str ());
      a_histo_gamma_energy_max.fill(E_max);
    }

  return dpp::base_module::PROCESS_OK;
}

bool snemo_gamma_tracking_efficiency_module::_compare_sequences(const gamma_dict_type & simulated_gammas_,
                                                                const gamma_dict_type & reconstructed_gammas_)
{
  _efficiency_.nevent++;

  if (reconstructed_gammas_.empty() && simulated_gammas_.empty())
    {
      DT_LOG_DEBUG(get_logging_priority(), "No gammas have been catched and reconstructed !");
      _efficiency_.nmiss++;
      return false;
    }

  _efficiency_.ntotal += simulated_gammas_.size();

//  std::cout << "simulated gamma size " <<simulated_gammas_.size() << std::endl;

  if(simulated_gammas_.size()>0)
    {
      _efficiency_.nevent_gammas++;
    }

  if (simulated_gammas_.size() > 1) {
    //*    DT_LOG_WARNING(datatools::logger::PRIO_WARNING, "More than one gammas simulated");
  }
  if (simulated_gammas_.size() == 0) {
    DT_LOG_WARNING(datatools::logger::PRIO_WARNING, "No gammas simulated");
  }

  if (get_logging_priority() >= datatools::logger::PRIO_DEBUG) {
    DT_LOG_DEBUG(get_logging_priority(), "Simulated gammas :");
    for (auto i : simulated_gammas_) {
      std::ostringstream oss;
      oss << "Gamma #" << i.first << " :";
      for (auto icalo : i.second) {
        oss << " -> " << icalo;
      }
      DT_LOG_DEBUG(get_logging_priority(), oss.str());
    }
    DT_LOG_DEBUG(get_logging_priority(), "Reconstructed gammas :");
    for (auto i : reconstructed_gammas_) {
      std::ostringstream oss;
      oss << "Gamma #" << i.first << " :";
      for (auto icalo : i.second) {
        oss << " -> " << icalo;
      }
      DT_LOG_DEBUG(get_logging_priority(), oss.str());
    }
  }

  size_t tmp_ngood_gammas = 0 ;
  for (auto irec : reconstructed_gammas_) {
    const calo_list_type & a_rec_list = irec.second;
    for (auto isim : simulated_gammas_) {
      const calo_list_type & a_sim_list = isim.second;
      if (a_rec_list.size() != a_sim_list.size()) continue;
      const bool are_same = std::equal(a_rec_list.begin(), a_rec_list.end(), a_sim_list.begin());
      if (are_same) {
        _efficiency_.ngood++;
        tmp_ngood_gammas++;
        DT_LOG_DEBUG(get_logging_priority(), "Sequences are identical !");
        break;
      }
      //* else
        //* DT_LOG_WARNING(get_logging_priority(), "Sequences are different !");
    }
  }

  if(tmp_ngood_gammas == simulated_gammas_.size() && simulated_gammas_.size() > 0)
    {
      DT_LOG_DEBUG(get_logging_priority(), std::endl << "°°°°°° Fully good event with at least one gamma ! °°°°°°" << std::endl);
      _efficiency_.ngood_event++;
      return true;
    }
  else
    {
      if(simulated_gammas_.size() > 0)
        DT_LOG_WARNING(get_logging_priority(), std::endl<<" ||||||||||||||||" << std::endl << "tmp_ngood_gammas " << tmp_ngood_gammas << std::endl
                       << "simulated_gammas_.size() " << simulated_gammas_.size() << std::endl
                       << "reconstructed_gammas_.size() " << reconstructed_gammas_.size() << std::endl );
      return false;
    }
}
bool snemo_gamma_tracking_efficiency_module::_compare_sequences_cluster(const gamma_dict_type & simulated_gammas_,
                                                                const gamma_dict_type & clustered_gammas_)
{
  if (clustered_gammas_.empty() && simulated_gammas_.empty())
    {
      DT_LOG_DEBUG(get_logging_priority(), "No gammas have been catched and clustered !");
      return false;
    }

  if(simulated_gammas_.size()>0)
    {
      _no_gt_efficiency_.no_gt_nevent_gammas++;
    }

  //* if (simulated_gammas_.size() > 1) {
  //   DT_LOG_WARNING(datatools::logger::PRIO_WARNING, "Cluster : More than one gammas simulated");
  // }

  if (simulated_gammas_.size() == 0) {
    DT_LOG_WARNING(datatools::logger::PRIO_WARNING, "Cluster : No gammas simulated");
  }

  if (get_logging_priority() >= datatools::logger::PRIO_DEBUG) {
    DT_LOG_DEBUG(get_logging_priority(), "Simulated gammas :");
    for (auto i : simulated_gammas_) {
      std::ostringstream oss;
      oss << "Gamma #" << i.first << " :";
      for (auto icalo : i.second) {
        oss << " -> " << icalo;
      }
      DT_LOG_DEBUG(get_logging_priority(), oss.str());
    }
    DT_LOG_DEBUG(get_logging_priority(), "Clustered gammas :");
    for (auto i : clustered_gammas_) {
      std::ostringstream oss;
      oss << "Gamma #" << i.first << " :";
      for (auto icalo : i.second) {
        oss << " -> " << icalo;
      }
      DT_LOG_DEBUG(get_logging_priority(), oss.str());
    }
  }

  size_t tmp_ngood_gammas = 0;
  for (auto irec : clustered_gammas_) {
    const calo_list_type & a_rec_list = irec.second;
    for (auto isim : simulated_gammas_) {
      const calo_list_type & a_sim_list = isim.second;
      if (a_rec_list.size() != a_sim_list.size()) continue;
      const bool are_same = std::equal(a_rec_list.begin(), a_rec_list.end(), a_sim_list.begin());
      if (are_same) {
        tmp_ngood_gammas++;
        DT_LOG_DEBUG(get_logging_priority(), "Sequences are identical !");
        break;
      }
      // else
      //   DT_LOG_WARNING(get_logging_priority(), "Sequences are different !");
    }
  }

  if(tmp_ngood_gammas == simulated_gammas_.size() && simulated_gammas_.size() > 0)
    {
      DT_LOG_DEBUG(get_logging_priority(), std::endl << "°°°°°° Fully good event with at least one gamma ! °°°°°°" << std::endl);
      _no_gt_efficiency_.no_gt_ngood_event++;
      return true;
    }
  else
    {
      if(simulated_gammas_.size() > 0)
        DT_LOG_WARNING(get_logging_priority(), std::endl<<" &&&&&&&&&&&&&&&&" << std::endl << "tmp_ngood_gammas " << tmp_ngood_gammas << std::endl
                       << "simulated_gammas_.size() " << simulated_gammas_.size() << std::endl
                       << "clustered_gammas_.size() " << clustered_gammas_.size() << std::endl );
      return false;
    }
}

} // namespace analysis

  // end of snemo_gamma_tracking_efficiency_module.cc
  /*
  ** Local Variables: --
  ** mode: c++ --
  ** c-file-style: "gnu" --
  ** tab-width: 2 --
  ** End: --
  */
