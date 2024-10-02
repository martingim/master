% %choose which day to use the results from 
%%
lab_folder = "Lab_results/24_04_12/";
lab_surface_probes = ["f1425_a0_308.csv";
            "f1425_a0_308_r2.csv";
            "f1425_a0_308_r3.csv";
            "f1425_a0_616_r1.csv";
            "f1425_a0_154_r1.csv"];


%%
lab_folder = "Lab_results/24_09_18/";
lab_surface_probes = ["f1425_a308_r1.csv";
            "f1425_a308_r2.csv";
            "f1425_a308_r3.csv";
            "f1425_a616_r4.csv";
            "f1425_a616_r5.csv";
            "f1425_a616_r6.csv";
            "f1425_a154_r7.csv";
            "f1425_a154_r8.csv";
            "f1425_a154_r9.csv"];







%% plot the piston movements
figure;
hold on
for lab_folder_number = 1:size(lab_surface_probes, 1)
    run_folder = sprintf("%s%d/" , lab_folder, lab_folder_number);
    load(append(run_folder, "fil1.dat"));
    fil = fil1;
    t = 0.008:0.01:size(fil1, 1)/100;
    plot(t, fil-fil(1,1), DisplayName=sprintf("%d", lab_folder_number))
    legend()
end

%% Plot surface probes

for sensor=1:4
    figure;
    title(sprintf("surface probe:%d", sensor))
    hold on;
    for run_number=1:size(lab_surface_probes, 1)
        probes = load(append(lab_folder, lab_surface_probes(run_number)));
        plot(probes(:,sensor+2), DisplayName=sprintf("%d", run_number));
    end
    legend;
end

%% plot
for run_number=1:size(lab_surface_probes, 1)
    run_folder = sprintf("%s%d/" , lab_folder, run_number);
    figure;
    hold on
    title(sprintf("run number %d", run_number));
    load(append(run_folder, "fil1.dat"));
    fil = fil1;
    t = 0.008:0.01:size(fil1, 1)/100;
    plot(t, fil-fil(1,1))
    
    load(append(run_folder, "fil2.dat"));
    fil = fil2;
    plot(t, fil-fil(1,1))
    
    load(append(run_folder, "fil3.dat"));
    fil = fil3;
    plot(t, fil-fil(1,1))
    
    load(append(run_folder, "padle_ut.dat"));
    
    plot(t(1:size(padle_ut, 1)), (padle_ut-padle_ut(1)))
    xlim([0, 10])
    legend('fil1', 'fil2', 'fil3', 'padle ut')

end