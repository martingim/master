% This function is for Stokes theory substitution of trigonometric
% identities 
function funSimple = trigIden(f, theta)
funSimple = simplify(subs(f, (cos(theta))^2,(1+cos(2*theta))/2));
funSimple = simplify(subs(funSimple, (sin(theta))^2,(1-cos(2*theta))/2));
funSimple = simplify(subs(funSimple, cos(theta)*sin(theta),sin(2*theta)/2));
funSimple = simplify(subs(funSimple, (sin(theta))^3,(3*sin(theta)-sin(3*theta))/4));
funSimple = simplify(subs(funSimple, (cos(theta))^3,(cos(3*theta)+3*cos(theta))/4));
funSimple = simplify(subs(funSimple, sin(theta)*cos(2*theta), (sin(3*theta)-sin(theta))/2));
funSimple = simplify(subs(funSimple, cos(theta)*sin(2*theta), (sin(3*theta)+sin(theta))/2));
funSimple = simplify(subs(funSimple, (cos(theta))^4, (3+4*cos(2*theta)+cos(4*theta))/8));
funSimple = simplify(subs(funSimple, (sin(theta))^4, (3-4*cos(2*theta)+cos(4*theta))/8));
funSimple = simplify(subs(funSimple, (cos(theta)*sin(theta))^2, (1-cos(4*theta))/8));
funSimple = simplify(subs(funSimple, sin(theta)*cos(4*theta), (sin(5*theta)-sin(theta))/2));
funSimple = simplify(subs(funSimple, cos(theta)*sin(4*theta), (sin(5*theta)+sin(theta))/2));
end