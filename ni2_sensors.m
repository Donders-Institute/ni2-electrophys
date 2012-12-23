function sens = ni2_sensors(type)

switch type
  case 'eeg'
    % create an eeg electrode array
    [chanpos, tri] = icosahedron162;
    chanpos        = chanpos*10;
    chanpos(chanpos(:,3)<0,:) = [];
    
    sens.chanpos = chanpos;
    sens.elecpos = chanpos;
    for k = 1:size(chanpos,1)
      sens.label{k,1}    = sprintf('eeg%02d', k);
      sens.chantype{k,1} = 'eeg';
    end
    sens = ft_datatype_sens(sens);
    
  case {'meg_mag' 'meg'}
    % create an meg magnetometer array
    [chanpos, tri] = icosahedron642;
    chanpos        = chanpos*12.5;
    chanpos(:,3)   = chanpos(:,3) - 1.5;
    
    [x, y, z, chanpos, tri] = intersect_plane(chanpos, tri, [0 0 0], [1 0 0], [0 1 0]);
    coilori                 = normals(chanpos, tri, 'vertex');
  
    nchan        = size(chanpos,1);
    sens.chanpos = chanpos;
    sens.chanori = coilori;
    sens.coilpos = chanpos;
    sens.coilori = coilori;
    sens.tra     = eye(nchan);
    for k = 1:nchan
      sens.label{k,1}    = sprintf('meg%03d', k);
      sens.chantype{k,1} = 'megmag';
    end
    sens = ft_datatype_sens(sens);
    sens.tri = tri;
    
  case 'meg_grad_axial'
    % create an meg axial gradiometer array
    [chanpos, tri] = icosahedron642;
    chanpos        = chanpos*12.5;
    coilori        = normals(chanpos, tri, 'vertex');
    coilori(chanpos(:,3)<1.5,:) = [];
    chanpos(chanpos(:,3)<1.5,:) = [];
    chanpos(:,3)   = chanpos(:,3) - 1.5;
    
    nchan = size(chanpos,1);
    triin = sum(tri<=nchan,2)==3;
    tri   = tri(triin,:);
    
    sens.chanpos = chanpos;
    sens.chanori = coilori;
    sens.coilpos = [chanpos; chanpos+2.5*coilori];
    sens.coilori = [coilori; -coilori];
    sens.tra     = [eye(nchan) eye(nchan)];
    for k = 1:nchan
      sens.label{k,1}    = sprintf('meg%03d', k);
      sens.chantype{k,1} = 'megmag';
    end
    sens = ft_datatype_sens(sens);
    sens.tri = tri;
    
    
  case 'meg_grad_planar'
  otherwise
    error('unsupported type %s', type);
end
 