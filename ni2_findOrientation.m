function source1orient = ni2_findOrientation(source3orient,grid)
% function source1orient = ni2_findOrientation(source3orient)
%
% source3orient:  3 x Source_inside_pos x Time
% grid            .inside and .outside fields
%
% source1orient:  Source_all_pos x Time

source1orient=nan(size(grid.pos,1),size(source3orient,3));
for ii=1:length(grid.inside)
  [uu,ss,vv]= svd(squeeze(source3orient(:,ii,:))*squeeze(source3orient(:,ii,:))');
  source1orient(grid.inside(ii),:) = uu(:,1)'*squeeze(source3orient(:,ii,:));
end
