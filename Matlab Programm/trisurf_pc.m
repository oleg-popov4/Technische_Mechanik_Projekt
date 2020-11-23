function hh = trisurf_pc(t, px, py, data)
 
 Px = px(t(:));
 Py = py(t(:));
 T = reshape(1:length(Px), [size(t,1) size(t,2)]);
 Data = repmat(data, size(t,2),1);
 h = trisurf(T, Px, Py, Data);
 if nargout == 1, hh = h; end