function [ gq ] = CreateGQScheme(N)
%CreateGQScheme Creates GQ Scheme of order N
%   Creates and initialises a data structure 
    gq.npts = N;
    if (N > 0) && (N < 6)
        %order of quadrature scheme i.e. %number of Gauss points
        gq.gsw = zeros(N,1); %array of Gauss weights
        gq.xipts = zeros(N,1); %array of Gauss points
        switch N
            case 1  
              gq.gsw(1) = 2; 
              gq.xipts(1) = 0;
            case 2 % basically just sets n = 2
              gq.gsw(1) = 1;  
              gq.gsw(2) = 1;  
              gq.xipts(1) = -sqrt(1/3);
              gq.xipts(2) = sqrt(1/3);
            case 3
              gq.gsw(1) = 5/9;
              gq.gsw(2) = 8/9;
              gq.gsw(3) = 5/9;

              gq.xipts(1) = -sqrt(3/5);
              gq.xipts(2) = 0;
              gq.xipts(3) = sqrt(3/5);
            case 4
              gq.gsw(1) = (18 + sqrt(30))/36;
              gq.gsw(2) = (18 + sqrt(30))/36;
              gq.gsw(3) = (18 - sqrt(30))/36;
              gq.gsw(4) = (18 - sqrt(30))/36;

              gq.xipts(1) = -sqrt(3/7 -(2/7)*sqrt(6/5));
              gq.xipts(2) = sqrt(3/7 -(2/7)*sqrt(6/5));
              gq.xipts(2) = -sqrt(3/7 +(2/7)*sqrt(6/5));
              gq.xipts(4) = sqrt(3/7 +(2/7)*sqrt(6/5));

            case 5
              gq.gsw(1) = 0;
              gq.gsw(2) = -(1/3)*sqrt(5 - 2*sqrt(10/7));
              gq.gsw(3) = (1/3)*sqrt(5 - 2*sqrt(10/7));
              gq.gsw(4) = -(1/3)*sqrt(5 + 2*sqrt(10/7));
              gq.gsw(5) = (1/3)*sqrt(5 + 2*sqrt(10/7));

              gq.xipts(1) = 128/225;
              gq.xipts(2) = (322 + 13*sqrt(70))/900;
              gq.xipts(2) = (322 + 13*sqrt(70))/900;
              gq.xipts(4) = (322 - 13*sqrt(70))/900; 
              gq.xipts(5) = (322 - 13*sqrt(70))/900;
        end
              
    else
      fprintf('Invalid number of Gauss points specified');
    end 
end

