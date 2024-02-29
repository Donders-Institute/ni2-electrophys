[activity, time] = ni2_activation;
figure; plot(time, activity);

%%

[activity, time] = ni2_activation('frequency', 100);
figure; plot(time, activity);

%%

[activity, time] = ni2_activation('frequency', 5);
figure; plot(time, activity);

%%

headmodel = ni2_headmodel('type', 'singleshell');
figure; ft_plot_headmodel(headmodel, 'axes', true);

%%

headmodel = ni2_headmodel('type', 'spherical', 'nshell', 3);
figure; ft_plot_headmodel(headmodel, 'axes', true);

%%

headmodel = ni2_headmodel('type', 'spherical', 'nshell', 1);
figure; ft_plot_headmodel(headmodel, 'axes', true);

%%

sens = ni2_sensors('type', 'eeg');
figure; ft_plot_sens(sens, 'label', 'number')

%%

sens = ni2_sensors('type', 'ctf151');
figure; ft_plot_sens(sens, 'label', 'label', 'coil', false)
figure; ft_plot_sens(sens, 'label', 'label', 'coil', true)

%%

sens = ni2_sensors('type', 'meg');
figure; ft_plot_sens(sens, 'label', 'number', 'coil', true)

%%

leadfield = ni2_leadfield(sens, headmodel, [0 0 6 1 0 0]);
figure; ni2_topoplot(sens, leadfield)

%%

data = leadfield * activity;

figure; ni2_topomovie(sens, data, time)

%%

sourcemodel = ni2_sourcemodel('type', 'mesh');

figure; ft_plot_mesh(sourcemodel)
figure; ft_plot_mesh(sourcemodel.pos(sourcemodel.inside,:))

%%

sourcemodel = ni2_sourcemodel('type', 'grid', 'resolution', 1);

figure; ft_plot_mesh(sourcemodel)
figure; ft_plot_mesh(sourcemodel.pos(sourcemodel.inside,:))
