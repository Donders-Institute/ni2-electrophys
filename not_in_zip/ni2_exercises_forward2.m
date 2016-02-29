% page 1
diploc = [0 0 8];
sens = ni2_sensors('type','eeg');
headmodel = ni2_headmodel('type', 'spherical', 'nshell', 3);
leadfield = ni2_leadfield(sens, headmodel, diploc);

[data, time] = ni2_activation('ncycle', 0.5, 'frequency', 1);

% page 2
dipmom = [1;0;0];
sensordata = leadfield*dipmom*data;

ni2_topomovie(sens, sensordata, time);

[data2, time] = ni2_activation('ncycle', 0.5, 'frequency', 1, 'latency', 0.6);
data3 = zeros(size(data2));
sourcedata = [data; data2; data3];
figure;plot(time, sourcedata);

sensordata2 = leadfield*sourcedata;
ni2_topomovie(sens, sensordata2, time);

% page 3
[data3, time] = ni2_activation('ncycle', 0.4, 'frequency', 1, 'latency', 0.6);
sourcedata = [data; data2; data3];
dippar = [0 0 8.5 0 1 0; 3 5 6 0.7 0 -0.7; -5 6 3 0.5 -0.5 0.6];
leadfield = ni2_leadfield(sens, headmodel, dippar);
sensordata = leadfield*sourcedata;
ni2_topomovie(sens, sensordata, time);

% page 4
dippar = [0 0 6 1 0 0];
headmodeleeg = ni2_headmodel('type', 'spherical', 'nshell', 3);
senseeg = ni2_sensors('type', 'eeg');
leadfieldeeg = ni2_leadfield(senseeg, headmodeleeg, dippar);
ni2_topoplot(senseeg, leadfieldeeg);

headmodelmeg = ni2_headmodel('type', 'spherical', 'nshell', 1);
sensmeg = ni2_sensors('type', 'meg');
leadfieldmeg = ni2_leadfield(sensmeg, headmodelmeg, dippar);
ni2_topoplot(sensmeg, leadfieldmeg);

% page 6
headmodel = ni2_headmodel('type', 'spherical', 'nshell', 3);
sens = ni2_sensors('type', 'eeg');
leadfield = ni2_leadfield(sens, headmodel, [0 0 6 1 0 0]);
[data, time] = ni2_activation;
sensordata = leadfield*data;

sensordataLM = sensordata - repmat(sensordata(74,:),[91 1]);

figure;ni2_topoplot(sens, sensordata(:,485));
figure;ni2_topoplot(sens, sensordataLM(:,485));

figure;plot(time, sensordata);
figure;plot(time, sensordataLM);

sensordataLMRM = sensordataLM - repmat(0.5.*sensordataLM(87,:),[91 1]);
