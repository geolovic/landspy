dem = GRIDobj('small25.tif');
fd = FLOWobj(dem);
% fac = flowacc(fd);
% w = fac > 1000;
% s = STREAMobj(fd, w);
% fac = fac.Z;
% save 'mlab_files/fac_tunez.mat' fac;
% fd.ixcix  = zeros(fd.size,'uint32');
% fd.ixcix(fd.ix) = uint32(1):uint32(numel(fd.ix));
