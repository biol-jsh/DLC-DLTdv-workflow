function out = changeOfCoordinateSys(A, u, v, w, tijk)
% changeOfCoordinateSys gets an input matrix, then, based on unitary vectors 
% from both old and new coordinates systems gets transorms the matrix
% the old coordinate system (ijk) is assumed to be the identity matrix or global coord
% so that OUT is based on uvk. tijk becomes the origin.
% ______________________________________________________________________________________ %
%
% INPUT: 
%   - A: input matrix of order [n, 3]
%   - u/v/w: vectors of order [n, 3] that stablish the new coordinate system
%            these must be orthogonal to each other and obey the right hand rule.
%   - tijk: translation vector of order [n,3] based on the ijk coordinate system
%           and centered at the new origin of uvw
%
% OUTPUT:
%   - out: coordinates in uvw frame of reference with tijk as the new origin
%
% FUTURE WORK:
%   - TODO: add input validation
%
%``````````````````````````````````````````````````````````````````````````````````````` %   
% -------------------------------------------------------------------------------------- %

  %get the size
  len = size(A,1);
  out = zeros(length(len),3);
  
  %assume base reference
  i = repmat([1,0,0], [len,1]);
  j = repmat([0,1,0], [len,1]);
  k = repmat([0,0,1], [len,1]);
  
  %ensure that they are unitary
  u = u./vecnorm(u,2,2);
  v = v./vecnorm(v,2,2);
  w = w./vecnorm(w,2,2);

  %first get A from the new origin but still in ijk
  tA = A - tijk;
  
  %transform every point
  for(ii=1:len)
    %make rotation matrix
    R = [dot(i(ii,:),u(ii,:)), dot(j(ii,:),u(ii,:)), dot(k(ii,:),u(ii,:)); ...
         dot(i(ii,:),v(ii,:)), dot(j(ii,:),v(ii,:)), dot(k(ii,:),v(ii,:)); ...
         dot(i(ii,:),w(ii,:)), dot(j(ii,:),w(ii,:)), dot(k(ii,:),w(ii,:));];
    
    %apply transformation
    out(ii,:) = R*(tA(ii,:)');
  end

end
% EOF EOF EOF EOF EOF EOF EOF EOF EOF EOF EOF EOF EOF EOF EOF EOF EOF EOF EOF EOF EOF EOF 