function funSimple = hyperIden(f, k, h, z)

% coefficients substitution,
funSimple = subs(f, (sinh(k*h))^2, (cosh(2*k*h)-1)/2);
funSimple = subs(funSimple, (cosh(k*h))^2, (cosh(2*k*h)+1)/2);
% coeff for \sinh^n
funSimple = subs(funSimple, (sinh(k*h))^3, sinh(k*h)*((cosh(2*k*h)-1)/2));
funSimple = subs(funSimple, (sinh(k*h))^4, ((cosh(2*k*h)-1)/2)^2);
funSimple = subs(funSimple, (sinh(k*h))^5, sinh(k*h)*((cosh(2*k*h)-1)/2)^2);
funSimple = subs(funSimple, (sinh(k*h))^6, ((cosh(2*k*h)-1)/2)^3);
funSimple = subs(funSimple, (sinh(k*h))^7, sinh(k*h)*((cosh(2*k*h)-1)/2)^3);
funSimple = subs(funSimple, (sinh(k*h))^8, ((cosh(2*k*h)-1)/2)^4);
funSimple = subs(funSimple, (sinh(k*h))^9, sinh(k*h)*((cosh(2*k*h)-1)/2)^4);
funSimple = subs(funSimple, (sinh(k*h))^10, ((cosh(2*k*h)-1)/2)^5);
funSimple = subs(funSimple, (sinh(k*h))^11, sinh(k*h)*((cosh(2*k*h)-1)/2)^5);
funSimple = subs(funSimple, (sinh(k*h))^12, ((cosh(2*k*h)-1)/2)^6);
funSimple = subs(funSimple, (sinh(k*h))^13, sinh(k*h)*((cosh(2*k*h)-1)/2)^6);
funSimple = subs(funSimple, (sinh(k*h))^14, ((cosh(2*k*h)-1)/2)^7);
% related with depth, 3rd order
funSimple = subs(funSimple, cosh(k*(5*h+3*z))+cosh(k*(h+3*z)), 2*cosh(3*k*(h+z))*cosh(2*k*h));
% related with depth, 4th order
funSimple = subs(funSimple, cosh(2*k*(5*h+2*z))+cosh(2*k*(h-2*z)), 2*cosh(6*k*h)*cosh(4*k*h+4*k*z));
funSimple = subs(funSimple, cosh(2*k*(3*h+2*z))+cosh(2*k*(h+2*z)), 2*cosh(2*k*h)*cosh(4*k*h+4*k*z));
funSimple =subs(funSimple, cosh(4*k*z)+cosh(4*k*(2*h+z)), 2*cosh(4*k*h)*cosh(4*k*h+4*k*z));
% related with depth, 5th order
funSimple = subs(funSimple, cosh(k*(5*h-3*z))+cosh(k*(11*h+3*z)), 2*cosh(8*k*h)*cosh(3*k*h+3*k*z));
funSimple = subs(funSimple, cosh(k*(7*h+3*z))+cosh(k*(1*h-3*z)), 2*cosh(4*k*h)*cosh(3*k*h+3*k*z));
funSimple = subs(funSimple, cosh(3*k*(3*h-z))+cosh(3*k*(5*h+z)), 2*cosh(12*k*h)*cosh(3*k*h+3*k*z));
funSimple = subs(funSimple, cosh(k*(7*h-3*z))+cosh(k*(13*h+3*z)), 2*cosh(10*k*h)*cosh(3*k*h+3*k*z));
funSimple = subs(funSimple, cosh(k*(11*h-3*z))+cosh(k*(17*h+3*z)), 2*cosh(14*k*h)*cosh(3*k*h+3*k*z));
funSimple = subs(funSimple, cosh(3*k*(h-z))+cosh(3*k*(3*h+z)), 2*cosh(6*k*h)*cosh(3*k*h+3*k*z));
funSimple = subs(funSimple, cosh(k*(3*h+5*z))+cosh(k*(7*h+5*z)), 2*cosh(2*k*h)*cosh(5*k*h+5*k*z));
funSimple = subs(funSimple, cosh(k*(1*h+5*z))+cosh(k*(9*h+5*z)), 2*cosh(4*k*h)*cosh(5*k*h+5*k*z));
funSimple = subs(funSimple, cosh(k*(11*h+5*z))+cosh(k*(1*h-5*z)), 2*cosh(6*k*h)*cosh(5*k*h+5*k*z));
funSimple = subs(funSimple, cosh(k*(3*h-5*z))+cosh(k*(13*h+5*z)), 2*cosh(8*k*h)*cosh(5*k*h+5*k*z));
funSimple = subs(funSimple, cosh(5*k*(h-z))+cosh(5*k*(3*h+z)), 2*cosh(10*k*h)*cosh(5*k*h+5*k*z));
end

