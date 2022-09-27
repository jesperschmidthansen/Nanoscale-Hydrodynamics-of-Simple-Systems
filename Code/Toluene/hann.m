
%%%%
%%usage wdata = hann(data)
%%
%% Multiplying a Hann window to a data set.
%% If data is a matrix the window is applied 
%% on each column. 
%%
%% Copyright (C), 2008, Jesper Hansen. 
%% See the COPYING file for license agreements.
%%%%

function wdata = hann(data)

  [nr nc] = size(data);

  if ( nr == 1 || nc == 1 )

    ldata = length(data);
    w = 0.5.*( 1 + cos(2*pi.*[1:1:ldata]./(2*ldata-1)) ); 

    wdata = zeros(ldata,1);
 
    %%Have to make this loop unfortunately
    for n=1:ldata
      wdata(n) = w(n).*data(n); 
    end
 
  else

    w = 0.5.*( 1 + cos(2*pi.*[1:1:nr]./(2*nr-1)) );
    wdata = zeros(nr, nc);
    for n=1:nc
      wdata(:,n) = w'.*data(:,n);
    end

  end
       
end
