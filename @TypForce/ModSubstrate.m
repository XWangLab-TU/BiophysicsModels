function [f] = ModSubstrate(f,s,varargin)

ip = inputParser;
ip.CaseSensitive = false;
ip.addRequired('f', @(x) isa(x,'TypForce'));
ip.addRequired('s', @(x) isa(x,'ModSubstrate'));
ip.addParameter('mex', [], @isobject);
ip.parse(f,s,varargin{:});
%----------------------------------------------------------------------------------------
