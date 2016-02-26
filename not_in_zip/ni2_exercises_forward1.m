% page 4
[data, time] = ni2_activation;
figure;plot(time, data);

[data2, time2] = ni2_activation('latency',0.4,'frequency',5);
figure;hold on;
plot(time,data);
plot(time2,data2,'r');

% page 5
datamix = data+data2;
figure;plot(time, datamix);

datacombined = [data;data2];
figure;plot(time, datacombined);

datamix1 = data+data2;
datamix2 = sum(datacombined,1);
datamix3 = [1 1]*datacombined;

mix = [1 1];
datamix3 = mix*datacombined;
assert(isequal(datamix1, datamix2));
assert(isequal(datamix1, datamix3));
assert(isequal(datamix2, datamix3));

% page 6
mix = [0 1;0.1 0.9;0.25 0.75;0.5 0.5;0.75 0.25;0.9 0.1;1 0];
datamix = mix * datacombined;
figure;plot(time, datamix+repmat((1:7)',[1 1000]));

% page 8
headmodel = ni2_headmodel('type', 'spherical', 'nshell', 3);
sens = ni2_sensors('type', 'eeg');
figure;hold on;
ft_plot_vol(headmodel,'edgecolor','none');
ft_plot_sens(sens);

dippar1 = [0 0 6 1 0 0];
leadfield1 = ni2_leadfield(sens, headmodel, dippar1);
ni2_topoplot(sens, leadfield1); colorbar;

dippar2 = [0 0 6 2 0 0];
leadfield2 = ni2_leadfield(sens, headmodel, dippar2);
ni2_topoplot(sens, leadfield2); colorbar;

% page 9
ni2_topoplot(sens, leadfield1*2);
ni2_topoplot(sens, leadfield2);

dippar3 = [0 0 5.9 1 0 0;0 0 6.1 1 0 0];
leadfield3 = ni2_leadfield(sens, headmodel, dippar3);
topo=leadfield3*[1;1];
ni2_topoplot(sens, topo);colorbar

