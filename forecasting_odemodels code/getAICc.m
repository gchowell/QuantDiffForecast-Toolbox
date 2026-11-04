function [AICc,part1,part2,numparams]=getAICc(method1,dist1,params_num,fixI0,fval,n)
  
numparams=params_num;
    
numparams=numparams+fixI0;
    

if method1==0 & dist1==0 % Normal distribution -- one parameter for variance

    numparams=numparams+1;

elseif method1==0 & dist1==2 % Neg. binomial error structure where var = factor1*mean

    numparams=numparams+1;

elseif method1==3 | method1==4 %Neg. Binomial requires one more parameter (alpha)

    numparams=numparams+1;

elseif method1==5

    numparams=numparams+2;  % Neg. Binomial requires 2 more parameters (alpha,d)

elseif method1==6  % Laplace distribution -- one parameter for scale

    numparams=numparams+1;

end


switch method1
    
    case 0
        
        AICc= n*log(fval/n) + n*log(2*pi) + n + 2*numparams + (2*numparams*(numparams+1))/(n-numparams-1);
        
        part1= n*log(fval/n);
        
        part2= n*log(2*pi) + n + 2*numparams + (2*numparams*(numparams+1))/(n-numparams-1);
    
    case 6

        AICc= 2*n*log(2*fval/n) + 2*n + 2*numparams + (2*numparams*(numparams+1))/(n-numparams-1);
        
        part1 = 2*n*log(2*fval/n);
        
        part2=  2*n + 2*numparams + (2*numparams*(numparams+1))/(n-numparams-1);


    otherwise
        
        AICc=-2*(-fval) + 2*numparams + (2*numparams*(numparams+1))/(n-numparams-1);
        
        part1=-2*(-fval);
        
        part2=2*numparams + (2*numparams*(numparams+1))/(n-numparams-1);

end

%method1

%fval

%n

%numparams
