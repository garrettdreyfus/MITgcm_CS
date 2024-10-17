function runShelfDepthCheck()
  %% Depth at this seed is 465
  experiment_parameters = struct;
  experiment_parameters.tcline_deltaz = 100;
  experiment_parameters.trough_depth = 0;
  experiment_parameters.rand_topo = true;
  experiment_parameters.monitor_freq = 6;
  rng_seed = 16;
  experiment_parameters.rng_seed = rng_seed;
  experiment_parameters.random_amplitude = 0;
  experiment_parameters.saltflux = true;
  experiment_parameters.rbcs_temp = true;
  experiment_parameters.cavity_depth = -300;
  experiment_parameters.hydro_mode = 'cold';
  experiment_parameters.initial_state = 'cold';
  experiment_parameters.yicefront = 150;
  experiment_parameters.shelf_depth = 601;
  experiment_parameters.tcline_atshelf_depth = 0;
  currentFolder = pwd;
  experiment_parameters.cavity_width = 150;
  %widths = [75 250];
  fronts = [75,100,150];
  for i = 1:3
    experiment_parameters.yicefront = fronts(i);
    experiment_parameters.saltflux_factor = 0.50;
    path_part1 = convertStringsToChars(strcat("experiments/icefront",int2str(rng_seed),"/"));
    path_part2 = convertStringsToChars(strcat("coldersf50front",int2str(fronts(i))));
    full_path = convertStringsToChars(strcat("../",path_part1,path_part2));
    newexp(path_part1,path_part2,experiment_parameters);
    cd(full_path);
    system('sh upload_to_cluster.sh')
    cd(currentFolder);
  end
  %for i = 1:3
    %experiment_parameters.saltflux_factor = 0.00;
    %experiment_parameters.cavity_width = widths(i);
    %path_part1 = convertStringsToChars(strcat("experiments/width",int2str(rng_seed),"/"));
    %path_part2 = convertStringsToChars(strcat("sf10unsalted",int2str(widths(i))));
    %full_path = convertStringsToChars(strcat("../",path_part1,path_part2));
    %newexp(path_part1,path_part2,experiment_parameters);
    %cd(full_path);
    %system('sh upload_to_cluster.sh')
    %cd(currentFolder);
  %end
end
