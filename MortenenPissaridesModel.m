% MP_root
function MP_root=MortenenPissaridesModel(theta)
global kappa beta q_theta P mu z b theta delta 
MP_root= kappa./(beta.*q_theta)-P.*((1-mu)*(z-b)-kappa.*mu.*theta +(1-delta).*kappa./q_theta);
end