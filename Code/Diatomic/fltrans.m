
%%
%% usage:  F = fltrans(x, f, w)
%%  
%% Fourier-Laplace transform (one-sided Fourier transform)
%%
%%         / x_end
%%        | 
%% F(w) = | f(x) exp(-iwx) dx
%%        |
%%       / 0
%%
%% Input variables:
%%  f      : Data to be transformed (vector)
%%  x      : Free variable (vector)
%%  w      : Frequency (vector) 
%%
%% Returns Fourier-Laplace transform, F
%%
%% Copyright (C), 2014, Jesper Hansen. 
%% See the COPYING file for license agreements.
%%

function F = fltrans(x, f, freq)
 
  if( isvector(f)==false || isvector(x) == false )
    error('sft: 1 and 2 argument must be vectors');
  elseif ( length(x) ~= length(f) )
    error('sft: Length of first  and second argument must be the same');
  end
  
  l = length(freq);
  I=i;
  for n=1:l
    g = f.*exp(-I*freq(n).*x);
    F(n) = trapz(x, g);
  end
  
end
