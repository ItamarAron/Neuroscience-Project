function [out] = randinterval(s1,s2,s3,llim,ulim)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here
out= (ulim - llim)*rand(s1,s2,s3)+llim*ones(s1,s2,s3);



end

