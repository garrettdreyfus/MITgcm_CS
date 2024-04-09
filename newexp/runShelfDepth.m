function runShelfDepthCheck()
  %% Depth at this seed is 465
  experiment_parameters = struct;
  experiment_parameters.tcline_deltaz = 100;
  experiment_parameters.trough_depth = 0;
  experiment_parameters.rand_topo = true;
  experiment_parameters.monitor_freq = 6;
  rng_seed = 16;
  experiment_parameters.rng_seed = rng_seed;
  experiment_parameters.random_amplitude = 200;
  experiment_parameters.saltflux = true;
  experiment_parameters.rbcs_temp = false;
  experiment_parameters.cavity_depth = -300;
  experiment_parameters.cavity_width = 150;
  experiment_parameters.hydro_mode = 'warm';
  experiment_parameters.yicefront = 150;
  currentFolder = pwd;
  shelf_depths = [ 500 600 700];
  for i = 1:3
    d = shelf_depths(i);
    experiment_parameters.shelf_depth = d;
    experiment_parameters.tcline_atshelf_depth = 0;
    path_part1 = convertStringsToChars(strcat("experiments/sfwarm",int2str(rng_seed),"/"));
    path_part2 = convertStringsToChars(strcat("d",int2str(d)));
    full_path = convertStringsToChars(strcat("../",path_part1,path_part2));
    newexp(path_part1,path_part2,experiment_parameters);
    cd(full_path);
    system('sh upload_to_cluster.sh')
    cd(currentFolder);
  end
end
