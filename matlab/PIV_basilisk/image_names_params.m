function [image_params] = image_names_params()
%the name of the folder for one day at the lab
day_folder = [  "MEK4600_G3_21_10_02";
                "MEK4600_G3_21_10_02";
                "MEK4600_G3_21_10_02";
                "MEK4600_G3_21_10_02";
                "MEK4600_G3_21_10_02";
                "MEK4600_G3_21_17_02";
                "MEK4600_G3_21_17_02";
                "MEK4600_G3_21_17_02";
                "MEK4600_G3_21_17_02";
                "MEK4600_G3_21_17_02";
                "MEK4600_G3_21_24_02";
                "MEK4600_G3_21_24_02";
                "MEK4600_G3_21_24_02";
                "MEK4600_G3_21_24_02";
                "MEK4600_G3_21_24_02"];


%the name of the folder for the run
run_folder = ["2021-02-10_run1_C001H001S0001";
                "2021-02-10_run2_C001H001S0001";
                "2021-02-10_run3_C001H001S0001";
                "2021-02-10_run4_C001H001S0001";
                "2021-02-10_run5_C001H001S0001";
                "2021-02-17_run6_C001H001S0001";
                "2021-02-17_run7_C001H001S0001";
                "2021-02-17_run8_C001H001S0001";
                "2021-02-17_run9_C001H001S0001";
                "2021-02-17_run10_C001H001S0001";
                "2021-02-24_run11_C001H001S0001";
                "2021-02-24_run12_C001H001S0001";
                "2021-02-24_run13_C001H001S0001";
                "2021-02-24_run14_C001H001S0001";
                "2021-02-24_run15_C001H001S0001"];

%The number for two consecutive images with the wave in the middle of the
%image. With three pairs of images for each run.     %run number
image_pairs = [ "016" "017" "088" "089" "159" "160"; %1
                "019" "020" "054" "055" "089" "090"; %2
                "042" "043" "069" "070" "096" "097"; %3
                "007" "008" "033" "034" "060" "061"; %4
                "013" "014" "047" "048" "080" "081"; %5
                "097" "098" "130" "131" "165" "166"; %6
                "028" "029" "062" "063" "096" "097"; %7
                "015" "016" "040" "041" "069" "070"; %8
                "017" "018" "044" "045" "070" "071"; %9
                "015" "016" "049" "050" "084" "085"; %10
                "027" "028" "062" "063" "096" "097"; %11
                "041" "042" "075" "076" "109" "110"; %12
                "036" "037" "071" "072" "105" "106"; %13
                "010" "011" "037" "038" "064" "065"; %14
                "005" "006" "032" "033" "059" "060";]; %15

%Fill image_name with the names of the images
image_name = string(zeros(size(run_folder, 2), 6));
for i=1:size(run_folder, 1)
    for j=1:size(image_pairs, 2)
        image_name(i,j) = day_folder(i,:) + "/" + run_folder(i) + "/" + run_folder(i) +  "000" + image_pairs(i,j) + ".bmp";
    end
end

%the frequencies for the runs
frequency = [1.75;1.75;2.25;2.25;1.75;
            1.75;1.75;2.25;2.25;1.75;
            1.75;1.75;1.75;2.25;2.25];
%the time between the frames
time_between_frames = [1/120;1/60;1/60;1/60;1/60;
                        1/60;1/60;1/60;1/60;1/60;
                        1/60;1/60;1/60;1/60;1/60;];


surface_start_stop = [  1000 5000;%can check first ones
                        1000 2000;
                        800 2000;
                        800 2000;
                        200 2000;
                        300 2000;
                        700 2000;
                        400 2000;
                        650 2000;
                        1000 2000;
                        1000 2000;
                        1000 2000;
                        1000 2000;
                        1000 2000;
                        1200 2000];

water_depth = [ 0.533;
                0.533;
                0.533;
                0.533;
                0.533;
                0.533;
                0.533;
                0.533;
                0.533;
                0.533;
                0.533;
                0.533;
                0.533;
                0.533;
                0.533;];
tform = coord_config();

image_params = containers.Map;
image_params('image_name') = image_name;
image_params('frequency') = frequency;
image_params('time_between_frames') = time_between_frames;
image_params('surface_start_stop') = surface_start_stop;
image_params('water_depth') = water_depth;
image_params('tform') = tform; 
image_params('day_folder') = day_folder;
end
