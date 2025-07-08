function [n] = numDigitInt(N,varargin)
        % numDigitInt computes the number of digit of the given integer
        % input: 
        % N - must be an integer
ip = inputParser;
ip.CaseSensitive = false;
ip.addRequired('N', @(x) isinteger(x));
ip.parse(N,varargin{:});
%--------------------------------------------------------------------------
n=floor(log10(abs(double(N))+1)) + 1;