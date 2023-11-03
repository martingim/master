function piv = initpass(Uspline,Vspline)
% INITPASS
% 
% piv = initpass(piv)
% 
% piv = initpass(Uspline,Vspline)

  
  tic;    
    if(nargin<2)
        piv.Vspline = Uspline.Vspline;
        piv.Uspline = Uspline.Uspline;
    else    
        piv.Uspline = Uspline;
        piv.Vspline = Vspline;
    end
    
    piv.passes = 1;  
    piv.pass{piv.passes} = 'init';
    piv.opt{piv.passes} = [];
  
    fprintf('PIV init pass %d: ',piv.passes);  
 toc
    
  
