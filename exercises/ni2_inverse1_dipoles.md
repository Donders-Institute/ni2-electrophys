# 1 Fitting EEG/MEG activity with dipole models

This document contains the MATLAB exercises that form part of the course “Neuroimaging II” relating to dipole modelling, or dipole fitting.

## 1.1 Introduction

These exercises treat the topic of dipole modelling. This technique is a so-called parametric inverse method, because we will estimate a set of parameters yielding a model of the observed data, which is aimed at explaining the data in an optimal way. The modelling assumption we make here is that the spatial distribution of measured potential or magnetic field at any moment in time is 'caused’ by a limited number of equivalent current dipoles (ECDs). In the past sessions you have learnt that we can build models of the potential difference or magnetic field distribution for a given set of dipole parameters (3 location, and 3 moment parameters). This is the so-called forward model. The idea of dipole modelling is as follows: given a set of parameters we build a forward model of the data, and compare this forward model with the measured data, by computing a measure of goodness-of-fit. Then, we manipulate the parameters a tiny bit, recompute the forward model, and see whether the goodness-of-fit is improved. If the goodness-of-fit is improved, then the current model of the data apparently is better than the previous one. We continue changing the parameters until we cannot further improve the goodness-of-fit.

# 2 Fitting a model to the observed data

The first thing we need to do is to determine which data matrix we are going to model. As you know by now, the measured data can be represented in a spatio-temporal matrix, where the rows represent the channels, and the columns represent the spatial topography of the observed activity for a given time point. We will start by modelling the data that is observed at a single time point. Then, we need to make an assumption about the number of dipoles (and thus the number of parameters) that we believe underly the observed data. Here, we start with assuming that the data can be modelled with just a single dipole.

After doing these exercises

-   You understand how the observed data can be fitted to a model.
-   You understand how to quantify the goodness-of-fit

We begin by generating some sensor-level data, by now this hopefully starts to look familiar:

    [data, time] = ni2_activation;
    sens = ni2_sensors('type', 'eeg');
    headmodel = ni2_headmodel('type', 'spherical', 'nshell', 3);
    leadfield = ni2_leadfield(sens, headmodel, [4.9 0 6.2 0 1 0]);
    sensordata = leadfield_data+randn(91, 1000)*1e-3;

We have added some noise to the data, in order to make it look more realistic. Now we select a single time slice from the data matrix, which will serve as our observed topography:

    topo_observed = sensordata(:, 500);
    figure; ni2_topoplot(sens, topo_observed);

What is typically done in dipole fitting is that the model of the observed data is created in a two-step procedure. First, given a set of location parameters, the forward model is created, yielding an Nchannels x (Ndipoles x 3) leadfield matrix, where each triplet of columns in the leadfield matrix represents the topography of a unit-amplitude dipole at a given location, with an orientation along one of the axes of the coordinate system. In a second step, the parameters that describe the dipole(s’) moment are estimated using linear estimation. This can be seen as a linear regression problem, where the leadfield matrix is the 'design matrix’, and the observed topography is the dependent variable.

In MATLAB this can be done in the following way. We assume a dipole at an arbitrary location and first model the leadfield.

    leadfield = ni2_leadfield(sens, headmodel, [5 5 4]);

Next, we are going to estimate the remainder of the parameters, i.e. the 3 dipole moment parameters, using linear estimation/least-squares regression. In MATLAB this can be easily done using the `\` or backslash operator.

    dipmom = leadfield\topo_observed;

> Revisit the MATLAB-exercises from the advanced MATLAB course (section 5: solving systems of linear equations) to refresh your understanding of this estimation step. If you haven’t followed this particular course, ask your fellow students for a copy of this document.

The modelled data can now easily be obtained by multiplying the leadfield with the moment parameters.

> Do this multiplication and visualize the topographical distribution using the `ni2_topoplot` function. Create a variable called `topo_modelled` representing this topographical distribution.

Now we can quantify the difference between the observed and modelled data, by summing the squared differences, and relating this number to the sum of the squared data values. The smaller this number, the better the fit.

    sumsq = sum((topo_observed-topo_modelled).^2)./sum(topo_observed.^2);

The value we have obtained for sumsq does not mean much, unless it is compared to the sumsq obtained for a fitted model with different parameters (i.e. a dipole at another location and or with different orientation).

> Model a dipole at location [-5 5 4], estimate the dipole moment given the observed data, and compute the relative sum of squared differences. Compare this to the value obtained above. Is the model fit better or worse?

# 3 Finding the optimal model

Now we know how to model the observed data using the leadfields created for one or more dipoles with a prespecified location, and we know how to quantify the goodness-of-fit between the modelled and observed data. Next, we need to consider the strategies that can be used to find the optimal model. Since the leadfields are non-linear functions of the parameters, there is no easy analytic solution to this problem. Therefore, the implicit strategy is to sample the parameter space, and to quantify for each setting of the parameters the goodness-of-fit. The parameter settings yielding the best model fit are selected. Typically, it does not make sense to just start placing dipoles at random locations and to hope that you will find the model with the best overall fit (out of all possible models. Particularly, when more than one dipole is assumed, the number of potential combinations of locations quickly explodes, and becomes unmanageable. However, in the single dipole case this seems to be a reasonable strategy. Under the assumption that the error landscape (i.e. 1 – goodness-of-fit expressed as a function of dipole location) is relatively smooth, one can sample the total set of possible source locations on a 3-dimensional grid of equally spaced dipole locations, and select the location that yields the smallest error. From this location, one could start a non-linear search to find a final solution.
After these exercises

-   You will understand the concept of the error function
-   You will understand the concept of a grid search

To create an error function, one simply needs to repeat multiple times the model fitting described in the previous section. One sensible strategy is to create a regular 3-dimensional grid of dipole positions that can be used to sample the parameter space. Here we do this with the following function call:

    sourcemodel = ni2_sourcemodel('type', 'grid', 'resolution', 1);

As a little aside, this creates a variable in MATLAB that is a so-called structure, which is a special type of variable that is convenient when working with data objects that have multiple attributes, i.e. features that belong together. These features are stored in so-called fields, which contain the actual data. In this example, if you type `sourcemodel` on the command line, you will see:

    sourcemodel =
            pos: [3610x3 double]
         inside: [1662x1 double]
        outside: [1948x1 double]
            dim: [19 19 10]

which means that the structure called `sourcemodel` has 4 fields, called `pos`, `inside`, `outside`, and `dim`. The most relevant field is the pos-field which contains a matrix with the positions defined in 3D-space. Note, that the number of positions is 3610, which equals the product of the elements in the dim-field. The dim-field defines the number of dipoles in each of the three directions of the 3D-box that defines the search space. The inside and outside fields are vectors which index which of the dipole positions are within the volume conduction model, and which ones are outside it. Note, that the total number of elements in the inside and outside fields matches the total number of positions. FieldTrip makes extensive use of MATLAB-structures, so it is good to have some basic understanding of this type of variable.

For illustration purposes we will first look at a subset of all positions in this 3D-grid, and select those positions which have a z-coordinate of 6. In other words, we are selecting those points that are on a plane that goes almost through the location where we actually simulated the activity (we are cheating a bit here, because usually we don’t know this of course).

    pos = sourcemodel.pos(sourcemodel.inside,:);
    sel = find(pos(:, 3)==6);
    pos = pos(sel,:);

Now, what we can do is repeat the steps in the previous section for each of these points, i.e. we will model for each of the positions a dipole that optimally explains the observed topography and compute a measure of goodness-of-fit. We store these goodness-of-fit measures in a vector, so that we later can determine which position gave the best fit. Doing the same time multiple times can be easily solved with a for-loop.

> If you don’t know already about the concept of a for-loop, read a bit about this in the MATLAB-documentation (type ‘help for’ or ‘doc for’ for this purpose).

    sumsq = zeros(size(pos, 1), 1);
    for k=1:size(pos, 1)
      leadfield = ni2_leadfield(sens, headmodel, pos(k,:));
      dipmom = leadfield\\topo_observed;
      topo_modelled = leadfield*dipmom;
      sumsq(k)=sum((topo_observed-topo_modelled).^2)./sum(topo_observed.^2);
    end

> Verify that what’s happening within the for-loop is the same as you did in the previous section.

Now we have a variable sumsq that is a vector, rather than a single number. We can visualize this variable using the `plot` function.

> Do this.

We can now look for the position that yields the lowest sum of squared difference values, i.e. the one with the best model fit. We could achieve this by zooming in into the figure and writing down the x-coordinate where the sumsq-variable seems to have the lowest value, but we can also use MATLAB’s min function:

    [m, ix] = min(sumsq);
    pos(ix,:)

Depending on the noise in the data, this position will be more or less close to the actual simulated location. One can improve on this estimate by fitting a model to more than one time point. The idea is that the influence of noise on the single time point observed topographies is attenuated when combining across time points. Assuming a fixed dipole location, this extended temporal model is quite straightforward, we will return to this later.

First, for convenience, we will repeat the steps above by using a set of pre-computed leadfields, that have been computed for the full 3D-grid. Since the computation of the leadfields is relatively expensive (time-consuming) in terms of computation, it makes sense to compute them once when they will be used multiple times. Further down in this tutorial we will explore the effect of model assumptions on the outcome, and do this using the grid search, so it pays off to compute the leadfields. Although this is not a really difficult exercise (just do the following: leadfield = ni2_leadfield(sens, headmodel, pos)) the leadfields are in the ni2-folder where the matlab-code for this part of the course is located. When you type load('leadfields’) you will obtain a variable called leadfield which contain per location a 3-column leadfield matrix.

> What do these 3-columns per location mean?

Now we can very quickly compute the goodness-of-fit for all locations in the 3D-grid:

    pos = sourcemodel.pos(sourcemodel.inside,:);
    load('leadfields');
    sumsq = ones(size(pos, 1), 1);
    for k=1:size(pos, 1)
      ik=(k-1)*3+(1:3);
      dipmom = leadfield(:, ik)\\topo_observed;
      topo_modelled = leadfield(:, ik)*dipmom;
      sumsq(k, 1)=sum((topo_observed(:)-topo_modelled(:)).^2)./sum(topo_observed(:).^2);
    end

> Check which location yields the lowest value for sumsq.

> Verify whether you understand the line `ik=(k-1)*3+(1:3);`. Hint: you can specify a value for k on the command line, execute the command, and then explore the content of ik.

# 4 Adding the time dimension

When using the topography from a single time point for dipole modelling, noise in the data can negatively affect the result. To this end it makes sense to combine data across time points. Under the assumption that the model is stationary in terms of dipole location (and possibly also in terms of orientation) pooling across time points can be efficient in improving the model fit.

After these exercises:

-   You will appreciate the fact that including temporal information can improve the dipole fit.

For this exercise we will the simulated data that was generated in section 2. If you are not working through these exercises in a single session, please go back to section 2, and re-create the variable sensordata.

Previously, we have used time point '500’ for the extraction of the observed topography.

> Why do you think we used this particular time point? Hint: consider the parameters that were used to generate the simulated time course.

Yet, in our very simple generative model (with generative model we mean the model that actually underlies the simulated data) there is only one dipole with a fixed location and orientation. In other words, and this is something that was highlighted in a previous session, there is essentially only a single topography present in the data, which only varies as a function of time. Thus, in principle, using time point '499’ should yield the same result as using time point '501’, or, for that matter, any other time point. Obviously, since the strength of the activity varies over time, some time points are more informative than others.

> Bonus question: which time point would yield a better result: 475 or 550? Argue why. Hint: consider the parameters used to generate the simulated time course.

Now we are going to fit a single dipole model to a spatial topography from a single time point, just as in the previous section.

    topo_observed = sensordata(:, 500);
    pos = sourcemodel.pos(sourcemodel.inside,:);
    load('leadfields');
    sumsq = ones(size(pos, 1), 1);
    for k=1:size(pos, 1)
      ik=(k-1)*3+(1:3);
      dipmom = leadfield(:, ik)\\topo_observed;
      topo_modelled = leadfield(:, ik)*dipmom;
      sumsq(k, 1)=sum((topo_observed(:)-topo_modelled(:)).^2)./sum(topo_observed(:).^2);
    end
    [m, ix] = min(sumsq);
    pos(ix,:)

Note the use of (:) in the line of code just before the end statement. This is the instruction for MATLAB to reshape any matricial variable into a vector. As such this is not needed in case the topography consists of just a single time point (because the observed topography is already a vector), but we will need it when fitting multiple time points.

# 5 The real deal

So far, for didactical purposes, we have constrained ourselves and evaluated the error-function based on a grid search. In this section we are going to use one of FieldTrip’s core functions ft_dipolefitting to do a proper dipole fit on our simulated EEG-data.

After this section

-   You understand that the non-linear search for optimal dipole parameters can yield a solution very close to the ground truth.

We start again by ensuring that we have the following variables in the MATLAB-workspace: sensordata, headmodel and sens. If you don’t have these variables available, please go back to section 2 and create these variables.
Next, we need to create some variables that are recognized by FieldTrip. Most FieldTrip functions work by providing 2 input arguments, a cfg-structure that contains parameters that determine the behaviour of the function, and one (or more) data-structures providing the data on which the function operates.

Let’s first create the data variable:

    data = [];
    data.avg = sensordata;
    data.time = time;
    data.label = sens.label;
    data.elec = sens;
    data.dimord = 'chan_time';
    And the cfg-structure:
    cfg = [];
    cfg.gridsearch = 'no';
    cfg.model = 'regional';
    cfg.vol = headmodel;
    cfg.latency = [0.49 0.51];
    cfg.nonlinear = 'yes';
    cfg.numdipoles=1;

Now we can call the function:

    dip = ft_dipolefitting(cfg, data);

> Read the documentation that belongs to `ft_dipolefitting` and try to understand the `cfg` you specified.

The output variable to `ft_dipolefitting` has a field `dip` containing information about the optimal model. In particular, have a look at `dip.dip.pos`, and `dip.dip.mom`. We can also visualize the modelled topography, and compare this to the observed topography. These are represented in `dip.Vmodel` and `dip.Vdata`, respectively.  

> Visualize the 11th column of the modelled and observed topographies using the `ni2_topoplot` function.

# 6 It’s all about the assumptions

One important issue in dipole modelling is that the prior assumptions critically constrain the final model (and thus the model fit). If these assumptions don’t coincide with what’s actually in the data, this can lead to erroneous interpretations. This can work in two directions, either the model is too simplistic (i.e. you assume too few dipoles), or too complicated (you assume too many dipoles).

After these exercises:

-   You understand that if your model assumptions violate the underlying data, the results are suboptimal.

First, we create some sensor-level data that contains two sources:

    [data, time] = ni2_activation;
    [data2, time] = ni2_activation('frequency', 11, 'latency', 0.48);
    sens = ni2_sensors('type', 'eeg');
    headmodel = ni2_headmodel('type', 'spherical', 'nshell', 3);
    leadfield = ni2_leadfield(sens, headmodel, [4.9 0 6.2 0 1 0]);
    leadfield2 = ni2_leadfield(sens, headmodel, [-5.3 0 5.9 1 0 0]);
    sensordata = leadfield_data+leadfield2_data2+randn(91, 1000)*.7e-3;

Now we will perform the same grid search as in section 4, using the 490th time point for the observed topography:

    topo_observed = sensordata(:, 490);
    sourcemodel = ni2_sourcemodel('type', 'grid', 'resolution', 1);
    pos = sourcemodel.pos(sourcemodel.inside,:);
    load('leadfields');
    sumsq = ones(size(pos, 1), 1);
    for k=1:size(pos, 1)
      ik=(k-1)*3+(1:3);
      dipmom = leadfield(:, ik)\\topo_observed;
      topo_modelled = leadfield(:, ik)*dipmom;
      sumsq(k, 1)=sum((topo_observed(:)-topo_modelled(:)).^2)./sum(topo_observed(:).^2);
    end

> Determine which dipole position yields the best goodness-of-fit. How does this position relate to the actual positions where activity was simulated?

Now we are going to use the FieldTrip function to fit a single dipole:

    data = [];
    data.avg = sensordata;
    data.time = time;
    data.label = sens.label;
    data.elec = sens;
    data.dimord = 'chan_time';

    cfg = [];
    cfg.gridsearch = 'no';
    cfg.model = 'regional';
    cfg.vol = headmodel;
    cfg.latency = [0.49 0.51];
    cfg.nonlinear = 'yes';
    cfg.numdipoles = 1;
    dip = ft_dipolefitting(cfg, data);

> Evaluate the position of the best fitting dipole, and compare this with the known locations of the sources.

We can also fit a model with two dipoles, this can be easily achieved by changing the `cfg.numdipoles` option into 2.

> Do this and evaluate the result.

As you may have noticed, the result is not particularly accurate. The reason for this is that the optimization algorithm got trapped in a local minimum of the error function. This is more likely to happen, the more complicated the underlying model (i.e. more free parameters lead to a high-dimensional error function with a complicated structure and potentially many local minima). We can however inform the fitting algorithm with dipole positions from which to start the non-linear search. If these starting positions are sufficiently close to the actual source positions, the algorithm will converge to the correct solution.

> Specify `cfg.dip.pos=[4 0 6;-4 0 6]` and verify that the result is now much better.
