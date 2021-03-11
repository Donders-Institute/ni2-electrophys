# 1 Modelling EEG/MEG activity using distributed sources

This document contains the MATLAB exercises that form part of the course “Neuroimaging II” relating to the minimum norm inverse methods for underdetermined systems.

## 1.1 Introduction

If we have measured data and a pre-computed leadfield (based on known sensor positions and a volume conductor model), then we would like to compute an estimate of the sources. Many methods use a linear matrix to compute this 'inverse solution’, as a linear weighting of the senors to compute activity at a source location.

Throughout this homework, the equation representing the data and model is `y = L*s+n`. Here, `y` is the sensor data, `L` is the known leadfield, `s` is the source amplitude and `n` is sensor noise.

# 2 Non-unique solutions for underdetermined systems: using the minimum norm of the source strength as additional constraint

After this section, you will

-   Understand that multiple solutions exist for an underdetermined system
-   Understand that the pseudo-inverse of the leadfield matrix gives the solution with the minimum-norm of the source power

## 2.1 Underdetermined linear systems of equations

First begin with a 'toy’ composite leadfield matrix. It represents 3 source locations and 2 sensors.

    L = [1 0 1;
         0 1 1];

Although it is an unrealistic example, it is useful to gain some intuition. Let’s assume the data we measured at one time point is

    y = [2 1]';

Now we want to solve for what the source s amplitude is, given this data `y` and given the known leadfield `L`. Assume no noise `n` for now, so the equation is `y = L*s`. This matrix multiplication essentially is a linear system of equations, where the number of equations is 2 (in general: the number of rows in the leadfield matrix), and the number of unknowns is 3 (in general: the number of columns in the leadfield matrix. This never has an unique solution. In other words, there are many solutions that satisfy this equation. Let’s first explicitly rewrite the matrix multiplication. Refresh your understanding of matrix multiplication if needed.

    y(1) = L(1, 1)*s(1) + L(1, 2)*s(2) + L(1, 3)*s(3); % equation 1
    y(2) = L(2, 1)*s(1) + L(2, 2)*s(2) + L(2, 3)*s(3); % equation 2

Filling in the numbers, we get:

    y(1) = 1*s(1) + 0*s(2) + 1*s(3) = 2;
    y(2) = 0*s(1) + 1*s(2) + 1*s(3) = 1;
    
> These are not valid lines of MATLAB code, you should check the equations in your head.

Note that equation 1 represents a matrix multiplication of the first row of the leadfield matrix with the source vector, and that equation 2 represents the matrix multiplication of the second row of the leadfield matrix with the source vector. Note, also, that we have more unknowns than equations. In other words, we could take an arbitrary value for, say, s(3), and we will still be able to find a valid solution to the linear system of equations. By analogy, you may remember from high school geometry, that equations with 3 unknowns describe a plane in 3-dimensional space, and that 2 planes (if they are not parallel) intersect in a line. This line represents the many valid solutions to the linear system of equations.

> Assume s(3) to be 0, and compute s(1) and s(2). Next, assume s(3) to be 3, and compute s(1) and s(2).

In this toy example, it is very straightforward to parametrize the solution. Let’s call the value that we take for s(3) the value `a`. Using this substition, and applying it to the equations above, we get:

    s(1) = 2–a;
    s(2) = 1-a;
    s(3) = a;

We will use this parametrization in the next section.

## 2.2 The 'best’ solution based on additional constraints: minimize norm.

In the example presented above, mathematically it is equally valid to take a value of 1 for a, as it is to take a value of 100. The former will yield source amplitudes of (1, 0, 1), whereas the latter will yield source amplitudes of (-98,-99, 100). It may be realistic to assume that the source activity measured at any given instant by MEG/EEG is primarily due to active sources that are moderately active (a=1). The alternative would be that all sources would be highly active (a=100). A mathematically convenient constraint is to minimize the total power of the sources across the brain (rather than taking the amplitude, we minimize for the amplitude squared). The total power of the sources across the brain is also called the L2-norm (the total amplitude of the sources across the brain is represented by the so-called L1-norm).

To compute the L2-norm, we take the square root of the sum of squares of each element in the vector of source amplitudes (vector over all source locations). Here’s an example of source amplitude over 3 source locations:

    s = [0 -1 2]';
    
> Verify that this is indeed yet another solution for the example in the previous section. To which value of a does this correspond?

The source power (L2 norm) can be computed in a variety of ways in Matlab:

    snorm1 = sqrt(s(1)^2+s(2)^2+s(3)^2)
    snorm2 = sqrt(sum(s.^2))
    snorm3 = sqrt(s'*s)

> Compute the source norm for the three source solutions (for values of `a` of 0, 3, and the solution to the previous exercise).  Which of the three has the minimum norm of source power?

Rather than trying out random values for `a` to find the solution that yields the lowest source norm, we can use the parametrization that we defined in the previous section:

    snorm = sqrt(sum(s.^2));

Filling in the parametrization we obtained in section 2.1, we get:

    snorm = sqrt((2-a)^2+(1-a)^2+a^2);

Rearranging the terms, we get:

    snorm = sqrt((4-4a+a^2)+(1-2a+a^2)+a^2);
    snorm = sqrt(3a^2-6a+5);

If we forget about the the square root, we need to minimize (3a^2-6a+5) which can be done in a variety of ways, e.g. using the [ABC-formula](https://en.wikipedia.org/wiki/Quadratic_formula) to find the minimum of a parabolic equation, or take the derivative of this equation, and find the value for `a`, where the derivative is 0.

## 2.3 Inverting leadfield matrices

For situations with more than just a few sensors and a few sources, it is too tedious to solve the system of equations for the unknown values by hand. Fortunately, we can use matrix algebra to solve large systems of equations, and let the computer do the work.

First, by mean of a detour, let’s have a look at the following equation:

    y = 3*s

If we want to solve this (very simple system of linear equations) for s (the source amplitude), assuming that we know y, we simply divide each side of the equation by 3, so we get:

    1/3 * y = s

Instead of dividing each side of the equation by 3, we can also say that we are multiplying each side of the equation by 1/3. 1/3 is also known as the 'inverse’ of 3, which can also be written as 3^-1. When we are dealing with matrices, the same logic applies:

    Y        =      L * S
    L^-1 * Y = L^-1 L * S 	
    L^-1 * Y =          S

Here `L^-1*L = I`, where the matrix `I` is the identity matrix, all zeros off the main diagonal and  all ones on the diagonal; this is the matrix equivalent of the number "1".

> Again, these are not valid lines of MATLAB code. Just follow the computations in your head.
 
The above equation however only works in a very limited number of cases. The reason is that a matrix inverse is only defined for square matrices (i.e. in our case we would need an equal number of channels and sources), where moreover the leadfield matrix must fulfill the mathematical property that it is actually invertible.

We can use the so-called pseudo-inverse of a matrix to get a solution to the linear system of equations. The term 'pseudo’ refers to the fact that the matrix you get after taking the pseudo-inverse behaves a bit like an inverse matrix, but not exactly.

The pseudo-inverse in this case can be computed in math formulation as:

    Lpinv = L' * (L*L')^-1

The MATLAB `pinv` function achieves the same:

    Lpinv = pinv(L)
    
> Using the leadfield L created in section 2.1, verify that these two methods give the same answer for Lpinv.
>
> Also try to invert L with the regular ‘inv’ command. What do you get?
>
> Investigate the outcome of Lpinv\*L and L\*Lpinv.

In contrast to a proper inverse matrix, where both `A*A^-1 = I` and `A^-1*A = I`, the pseudoinverse exists in only one direction. This can be easily seen when inspecting `L*Lpinv` and `Lpinv*L`. The first one does not result in an identity matrix. The second equation results in an identity, because it is equivalent to `(L*L')^-1 * (L*L')`. The matrix products between the brackets are the same, both form a square matrices and they happen to be in general invertible, thus the product between the inverse of that product with itself will yield `I`.

## 2.4 Pseudo-inverse of leadfield gives minimum norm

The pseudo-inverse is not only a clever mathematical way of computing a new matrix that, when multiplied with the original in the correct order, gives back the identity matrix. It happens to be a solution to the undertermined linear system of equations that yields the minimum norm.

This can be proved as follows (note that we start assuming that we have the solution that yields the minimum norm):

We know that there are many possible solutions for the equation `y = L*s`. Let us use a capital `S` from now on to facilitate formatting; we will use `Sx` to indicate s<sub>x</sub> as `s` with the subscript `x`. If we just take two possible solutions, `Smn` and `Sx`, where `Smn` happens to be `L'*(L*L')^-1 * y` and `Sx` being a random other solution, we know that both `L*Smn` and `L*Sx` yield `y`. Thus, `L*Smn = L*Sx`, or equivalently `L*(Smn-Sx) = 0`.

Remember from section 2.2 that the norm of the solution can be computed from `s'*s`. This latter equation computes for each of the sources the square of the amplitude, and sums across sources (to get the norm, we have to take the square root of the result, but let’s not do this for now, and look at the squared norm). Likewise, we can do `(Sx-Smn)'*Smn`, where we compute for each of the sources the product between the amplitude modelled in Smn and the amplitude difference between s and Smn, and sum this across sources. Substituting `L'*(L*L')^-1 * y` for the second `Smn` in the equation, we get:

      (Sx-Smn)' * Smn                  =
      (Sx-Smn)' * L'   * (L*L')^-1 * y =
    ( (Sx-Smn)' * L')  * (L*L')^-1 * y =
    (  L * (Sx-Smn) )' * (L*L')^-1 * y =		
         0             * (L*L')^-1 * y =	0

Shuffling the brackets, we focus on the first 2 terms, where we can use the general matrix property that `A'B'=(BA)'`.

Above we already concluded that `L*(Sx-Smn)` is 0, so we can conclude that  `(Sx-Smn)' * Smn` is 0. This is because we can fill in a 0 in the last equation, and 0 times something else will be 0.

Now considering the squared norm of `Sx`, which is `Sx'*Sx`, we can apply a little trick: instead of using `Sx`, we use `((Sx-Smn)+Smn)`. The latter is of course exactly the same as `Sx`.

      Sx'           *   Sx                                          =
    ((Sx-Smn)+Smn)' * ((Sx-Smn)+Smn)                                =
     (Sx-Smn)'*(Sx-Smn) + (Sx-Smn)'*Smn + Smn'*(Sx-Smn) + Smn'*Smn

Using our previous result `(Sx-Smn)'*Smn=0`, we get

    (Sx-Smn)'*(Sx-Smn) + (s-Smn)'*Smn + Smn'*(Sx-Smn) + Smn'*Smn    =
    (Sx-Smn)'*(Sx-Smn) + 0            + 0             + Smn'*Smn    =
    (Sx-Smn)'*(Sx-Smn) +                                Smn'*Smn

or

    Sx'*Sx = Smn'*Smn + (Sx-Smn)'*(Sx-Smn)
    
The last equation tells us that the squared norm of `Sx` is always the sum of the squared norm of `Smn` and the squared norm of the difference between `Sx` and `Smn`.

Since squared numbers are always positive (or 0), we can conclude that the norm of `Sx` is always larger or equal to the norm of `Smn`. Thus, `Smn` represents the solution with the minimum norm.

# 3 The real deal

After this section you will

-   Have a basic understanding of how the minimum norm reconstruction works in practice.
-   Understand why the minimum norm reconstruction tends to overemphasize activity from superficial sources.
-   Know how to counteract the tendency for overemphasizing the superficial sources.
-   Understand that noise in the data projects onto the estimated sources.
-   Know how to counteract the contamination of the estimated source activity by the noise.

## 3.1 Start simple: simulated data without noise

We start by simulating some MEG data that contains two active sources.

    [data1, time1] = ni2_activation;
    [data2, time2] = ni2_activation('frequency', 11, 'latency', 0.48);
    sens = ni2_sensors('type', 'meg');
    headmodel = ni2_headmodel('type', 'spherical', 'nshell', 1);
    leadfield1 = ni2_leadfield(sens, headmodel, [4.9 0 6.2 0 1 0]); % position 2352 in grid
    leadfield2 = ni2_leadfield(sens, headmodel, [-5.3 0 5.9 1 0 0]); % position 2342 in grid
    sensordata = leadfield1_data1+leadfield2_data2;

Try and understand the steps above. Pay particular attention to the parameters of the simulated dipoles.

We now proceed to generate a MATLAB data-structure that FieldTrip understands. This data-structure is a collection of MATLAB-variables, organized in so-called fields, that belong together. An important aspect of these FieldTrip data structures is that the numeric data that is represented (in our case in the 'avg’ field) is accompanied by all information necessary to interpret this numeric data. For example, there is a field called 'time’, that indicates each time sample in seconds (i.e. it maps the columns of the 'avg’ field on a physical time axis). The 'label’ field specifies the name of each channel (and tells us which row in the 'avg’ field belongs to which channel).

    data        = [];
    data.avg    = sensordata;
    data.time   = time1;
    data.label  = sens.label;
    data.grad   = sens;
    data.cov    = eye(numel(sens.label));
    data.dimord = 'chan_time';

Next we will make a source reconstruction using the 'mne’ method of FieldTrip’s ft_sourceanalysis function. Before we can do this, we need to define our source model, i.e. the set of locations that we assume to be active. For now we assume that the active dipoles are distributed on a 3D regular grid, with a spacing of 1 cm between the dipoles:

    sourcemodel = ni2_sourcemodel('type', 'grid', 'resolution', 1);

    cfg                    = [];
    cfg.grid               = sourcemodel;
    cfg.vol                = headmodel;
    cfg.method             = 'mne';
    cfg.mne.prewhiten      = 'yes';
    cfg.mne.scalesourcecov = 'yes';
    cfg.mne.lambda         = 0;
    cfg.keepleadfield      = 'yes';
    source = ft_sourceanalysis(cfg, data);

Let’s now have a look at the reconstructed source activity.

> First verify that  dipole 2352 and 2342 correspond to the locations that are closest to the locations at which actual activity was simulated. Hint: check the content of the rows in source.pos.  

For each dipole location in the distributed source model, the estimated activity is represented in the source.avg.mom field. We can easily use the MATLAB plot command to visualize this:

    figure; plot(source.time, source.avg.mom{2352}); legend({'x' 'y' 'z'});

> Why are there three lines in the figure? Does the dipole moment correspond with the simulated dipole moment?
>
> Also visualize dipole 2342. What do you see?
>
> Also visualize dipole 2347 (and inspect to what position this corresponds). What do you see?
>
> Visualize dipole 2713. What do you see? What is the position of this dipole and how does the reconstructed amplitude at this location related to the reconstructed amplitude at location 2352?

This simple noise-less example illustrates two important things. First, activity is 'smeared’ out over various dipole locations and orientations. Second, the estimated activity at a location closer to the sensors than the location at which activity was simulated has a higher amplitude than the activity estimated at the location where activity was simulated. We will return to this feature of the minimum norm reconstruction in a later section.

## 3.2 Simulated data with noise

Let’s now simulate MEG sensor data with added noise:

    [data1, time1] = ni2_activation;
    [data2, time2] = ni2_activation('frequency', 11, 'latency', 0.48);
    sens = ni2_sensors('type', 'meg');
    headmodel = ni2_headmodel('type', 'spherical', 'nshell', 1);
    leadfield1 = ni2_leadfield(sens, headmodel, [4.9 0 6.2 0 1 0]); % position 2352 in grid
    leadfield2 = ni2_leadfield(sens, headmodel, [-5.3 0 5.9 1 0 0]); % position 2342 in grid
    sensordata = leadfield1_data1+leadfield2_data2+randn(301, 1000)*.7e-10;

Create a FieldTrip data structure:

    data        = [];
    data.avg    = sensordata;
    data.time   = time1;
    data.label  = sens.label;
    data.grad   = sens;
    data.cov    = cov(randn(301, 1000)'*.7e-10);
    data.dimord = 'chan_time';

In the field 'cov’ we create a covariance matrix that was designed to represent the covariance of the noise in the data. This will be a relevant item when we will discuss noise regularisation.

    sourcemodel = ni2_sourcemodel('type', 'grid', 'resolution', 1);

Do the source reconstruction:

    cfg                    = [];
    cfg.grid               = sourcemodel;
    cfg.vol                = headmodel;
    cfg.method             = 'mne';
    cfg.mne.prewhiten      = 'yes';
    cfg.mne.scalesourcecov = 'yes';
    cfg.mne.lambda         = 0;
    cfg.keepleadfield      = 'yes';
    source_noise = ft_sourceanalysis(cfg, data);

> Explore the source activation time courses of the same dipoles as in the previous step (2352, 2342 etc). Describe what you observe.    

The `cfg.mne.lambda` option was set to 0. This means that the inverse solution fits all data perfectly, where the data not only includes the activity from the sources of interest, but also contains noise. We can use a non-zero lambda to compute a regularized minimum norm estimate, where this lambda parameter is used to quantify the contribution of the noise to the measured signals. The larger the value for lambda, the stronger the assumed noise. In combination with the regularization parameter, a regularized minimum-norm estimate also requires an estimate of the noise covariance matrix. This matrix represents the spatial structure in the noise. The noise covariance matrix can be estimated from the data, but experimenters sometimes also use an identity matrix. The latter strategy assumes implicitly that each channel in the data gets the same amount of uncorrelated noise.

Do the source reconstruction with regularisation:

    cfg                    = [];
    cfg.grid               = sourcemodel;
    cfg.vol                = headmodel;
    cfg.method             = 'mne';
    cfg.mne.prewhiten      = 'yes';
    cfg.mne.scalesourcecov = 'yes';
    cfg.mne.lambda         = 0.5;
    cfg.keepleadfield      = 'yes';
    source_noise_reg = ft_sourceanalysis(cfg, data);

> Explore the source activation time courses of the same dipoles as in the previous step (2351, 2342 etc). Describe what you observe.

Another way to explore the effect of regularisation is to look at the residuals of the model. This can be obtained in the following way:

> Repeat the regularisation, but now with a value for lambda of 0.1. Explore the source activation time courses of the same dipoles as in the previous step (2351, 2342 etc). Describe what you observe. Also apply a regularistion of 2.
>
> Describe the effects of an increasing regularisation parameter.


    L = cat(2, source_noise_reg.leadfield{source_noise_reg.inside});
    S = cat(1, source_noise_reg.avg.mom{source_noise_reg.inside});
    model = L*S;
    residual = sensordata-model;

> Compute and visualize the residuals for various values of lambda and describe what you see.

## 3.3 Minimum-norm estimates 'overestimate’ the amplitude of superficial sources

As we have seen in the previous sections, the minimum-norm estimate has a tendency to over-estimate the amplitude of dipoles that are close to the surface. This feature is a direct consequence of the minimum-norm constraint. In order to explain all measured data with a source model that has the lowest possible norm, the deep sources will be penalized because these need to have a strong activation in order to be picked up by the sensors in the first place.

    cfg                    = [];
    cfg.grid               = sourcemodel;
    cfg.vol                = headmodel;
    cfg.method             = 'mne';
    cfg.mne.prewhiten      = 'yes';
    cfg.mne.scalesourcecov = 'yes';
    cfg.mne.lambda         = 0.5;
    cfg.keepleadfield      = 'yes';
    cfg.normalize          = 'yes';
    cfg.normalizeparam     = 1;
    source_noise_lfnorm    = ft_sourceanalysis(cfg, data);
