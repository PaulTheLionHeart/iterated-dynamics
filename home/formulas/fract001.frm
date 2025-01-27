a---mand { ; Eli Brandt, (c) 1992
   z = c = pixel:
   z = z^c + sin(c);
   |z|<p1
}

a--mand { ; Eli Brandt, (c) 1992
   z = c = pixel:
   z = (sqr(z)+c) / |c|;
   |z|<p1
}

aa-c-to-z(XAXIS)  { ; Eli Brandt, (c) 1992
   z = pixel, c = z:
   z = c^z;
   |z| <= p1
}

aa-diamand { ; Eli Brandt, (c) 1992
   z = c = pixel:
   z = sqr(z)+c;
   (real(z)+imag(z))<p1
}

aa-mand-im(XAXIS) { ; Eli Brandt, (c) 1992
   z = c = pixel:
   z = sqr(z)+c;
   imag(z)<p1
}

aa-mand-re(XAXIS) { ; Eli Brandt, (c) 1992
   z = c = pixel:
   z = sqr(z)+c;
   real(z)<p1
}

aa-z-to-c(XAXIS) { ; Eli Brandt, (c) 1992
   z = pixel, c = z:
   z = exp(c*log(z));
   |z| <= p1
}

aa-z-to-z(XAXIS) { ; Eli Brandt, (c) 1992
   z = pixel, c = z:
   z = exp(z*log(z));
   |z| <= p1
}

aaa-deals { ; from DAN.FRM
   z = pixel, a=real(z), b=imag(z), i=((0-1)^0.5):
   x = real(z), y = imag(z), u = x*x + y*y + a,
   v = (0-2)*x*y + b,
   z = u + v*i;
   |z|<p1
}

aaa-dealsFix1 { ; from DAN.FRM
   ; fix by Ron Barnett
   z = pixel, a=real(z), b=imag(z), i=((0,-1)^0.5):
   x = real(z), y = imag(z), u = x*x + y*y + a,
   v = (0,-2)*x*y + b,
   z = u + v*i;
   |z|<p1
}

aaa-dealsFix2 { ; from DAN.FRM
   ; fix by Ron Barnett
   z = pixel, a=real(z), b=imag(z), i=((0,1)^0.5):
   x = real(z), y = imag(z), u = x*x + y*y + a,
   v = (0,2)*x*y + b,
   z = u + v*i;
   |z|<p1
}

alt (xaxis) {
   z=0, c=pixel, k=1:
   z=sqr(z) + c,
   c=c+k*p1/z, k=((11-3*k)*k-4)/2,
   |z| <= 4
}

cardioid {
   z=0, x=real(pixel), y=imag(pixel),
   c=x*(cos(y)+x*sin(y)):
   z=sqr(z)+c,
   |z| < 4
}

CGNewtonSinExp (XAXIS) { ; Chris Green
   ; Use floating point.
   z=pixel:
   z1=exp(z),
   z2=sin(z)+z1-z,
   z=z-p1*z2/(cos(z)+z1),
   .0001 < |z2|
}

ConformalMapping {
   c = pixel, RealZ = Real(c), ImagZ = Imag(c):
   RealZ = Sqr(RealZ) + (RealZ * ImagZ) + Real(c);
   ImagZ = Sqr(ImagZ) + (RealZ * ImagZ) + Imag(c);
   z = RealZ + (ImagZ * (0, 1)),
   |z| < 4
}

ConjMandelbrot(XAXIS) { ; Paul J. Horn
   ; This uses the last square function and should be the same as
   ; MandelConj, but with inside=bof60 and outside=summ and zooms
   ; it does not.
   z = Pixel, z = Sqr(conj(z)):
   z = z + Pixel
   z = Sqr(conj(z))
   LastSqr <= 4
}

CosInvZ(XYAXIS) {
   z=pixel,inv=1/pixel+p1:
   z=cos(inv/z),
   |z|<=4
}

CoshInvZ(XYAXIS) {
   z=pixel,inv=1/pixel+p1:
   z=cosh(inv/z),
   |z|<=4
}

{  ; Try p1=0, p2=4, fn1=sqr, fn2=exp, fn3=cosxx, for old DeepSpaceProbeTwo
   ; Try p1=0, p2=4, fn1=sqr, fn2=exp, fn3=log, for old Moth type
   ; Try p1=0, p2=4, fn1=sqr, fn2=cosxx, fn3=sin, for old ManInTheOzone type }

DeepSpaceProbe(XAXIS_NOPARM) { ; Jonathan Osuch
   ; Generalized by Tobey J. E. Reed
   ; Try p1=0, p2=4, fn1=sqr, fn2=sin, fn3=cosxx
   ; Note:  use floating point
   z   = p1, x   = 1:
   (x  <  10) * (z=fn1(z) + pixel),
   (10 <= x)  * (x<20)    * (z=fn2(z)+pixel),
   (20 <= x)  * (z=fn3(z) + pixel),
   x   = x+1,
   |z| <= p2
}

{  ; Try p1=0, p2=4, fn1=sqr, fn2=exp, fn3=cosxx, for old DeepSpaceProbeTwoC
   ; Try p1=0, p2=4, fn1=sqr, fn2=exp, fn3=log, for old MothC type
   ; Try p1=0, p2=4, fn1=sqr, fn2=cosxx, fn3=sin, for old ManInTheOzoneC type }

DeepSpaceProbeC(XAXIS_NOPARM) { ; Jonathan Osuch
   ; Generalized by Tobey J. E. Reed
   ; Try p1=0, p2=4, fn1=sin, fn2=exp, fn3=cosxx
   ; Note:  use floating point
   z    = p1, x    = 1:
   (z=fn1(z)+pixel)*(x<10)+(z=fn2(z)+pixel)*(10<=x)*(x<20)+(z=fn3(z)+pixel)*(20<=x),
   x   = x+1,
   |z| <= p2
}

dots { ; Eli Brandt, (c) 1992
   z = c = pixel:
   z = z^c + c,
   |z|<p1
}

DrChaosbrot { ; Michael Theroux
   ; fix and generalization by Ron Barnett
   ; more phi
   ; try p1 = 2.236067977 for the golden mean
   z = c = pixel:
   z = sqr(z) + (((p1 + 1)/2)+c)
   |z| <= 4
}

Element(xyaxis) { ; Michael Theroux
   ; fix and generalization by Ron Barnett
   ; phi lingam
   ; try p1 = 2.236067977 for the golden mean
   z = pixel:
   z = z*z*z*z + ((p1 + 1)/2)
   |z| <= 4
}

EvilEyeC(XAXIS_NOPARM) { ; Jonathan Osuch
   ; Generalized by Tobey J. E. Reed
   ; Try p1=1, p2=4, fn1=sin, fn2=cosxx
   ; Note:  use floating point
   z  =  p1,
   x  = |z|:
   (z  = fn1(z)+pixel) * (1<x)  +  (z=fn2(z) + pixel)*(x<=1),
   x  = |z|,
   x <= p2
}

Exipi (XAXIS) { ; Lee Skinner
   s = log(-1.,0.) / (0.,1.), z = Pixel:
   z = z ^ s + pixel,
   |z| <= 100
}

F'Liar1A { ; Generalization by Jon Horner of Chuck Ebbert formula.
   ; X: X is as true as Y
   ; Y: Y is as true as X is false
   ; Calculate new x and y values simultaneously.
   ; y(n+1)=abs((1-x(n) )-y(n) ), x(n+1)=1-abs(y(n)-x(n) )
   z = fn1(pixel):
   z = 1 - abs(imag(z)-real(z) ) + flip(1 - abs(1-real(z)-imag(z) ) ),
   |z| <= p1
}

F'Liar1C { ; Generalization by Jon Horner of Chuck Ebbert formula.
   ; X: X is as true as Y
   ; Y: Y is as true as X is false
   ; Calculate new x and y values simultaneously.
   ; y(n+1)=abs((1-x(n) )-y(n) ), x(n+1)=1-abs(y(n)-x(n) )
   z = fn1(pixel):
   z = 1 - abs(imag(z)-real(z) ) + flip(1 - abs(1-real(z)-imag(z) ) ),
   fn2(abs(z))<p1
}

F'Liar1D { ; Generalization by Jon Horner of Chuck Ebbert formula.
   ; X: X is as true as Y
   ; Y: Y is as true as X is false
   ; Calculate new x and y values simultaneously.
   ; y(n+1)=abs((1-x(n) )-y(n) ), x(n+1)=1-abs(y(n)-x(n) )
   z = fn1(pixel):
   z = p1 - abs(imag(z)-real(z) ) + flip(p2 - abs(1-real(z)-imag(z) ) ),
   |z| <1
}

F'M-SetInNewtonB(XAXIS) { ; use float=yes, periodicity=no
   ; set p1 >= 3, 1e-30 < p2 < .01
   z=0, c=fn1(pixel), cm1=c-1, cm1x2=cm1*2, twoop1=2/p1, p1xc=c*real(p1):
   oldz = z,
   z= (p1xc - z*cm1x2 )/( (sqr(z)*3 + cm1 ) * real(p1) ) + z*real(twoop1),
   |z - oldz| >= p2
}

flip0_man_j(ORIGIN) { ; Richard Hughes (Brainy Smurf)
   z=pixel:
   z = flip(sqr(z) + p1),
   |z| <= 4
}

flip0_man_m(XAXIS) { ; Richard Hughes (Brainy Smurf)
   z=0:
   z = flip(sqr(z) + pixel),
   |z| <= 4
}

flip1_man_j(ORIGIN) { ; Richard Hughes (Brainy Smurf)
   z=pixel, q = p1:
   q = flip(q),
   z = sqr(z) + q,
   |z| <= 4
}

flip1_man_m { ; Richard Hughes (Brainy Smurf)
   z=0, q = pixel:
   q = flip(q),
   z = sqr(z) + q,
   |z| <= 4
}

flip2_man_j(ORIGIN) { ; Richard Hughes (Brainy Smurf)
   z=pixel, q = p1:
   q = flip(q),
   z = flip(sqr(z) + q),
   |z| <= 4
}

flip2_man_m { ; Richard Hughes (Brainy Smurf)
   z=0, q = pixel:
   q = flip(q),
   z = flip(sqr(z) + q),
   |z| <= 4
}

flip3_man_j { ; Richard Hughes (Brainy Smurf)
   z = pixel:
   z = 1/flip(sqr(z) + p1),
   |z| <= 4
}

flip3_man_m(XAXIS) { ; Richard Hughes (Brainy Smurf)
   z = 0:
   z = 1/flip(sqr(z) + pixel),
   |z| <= 4
}

FlipLambdaM { ; Ron Barnett
   ; provides a "map" of locations for FlipLambdaJ
   ; Try "center-mag" with center = (0.49,0.31)
   ; mag = 10.4
   z = 0.5:
   z = pixel*z*(1-flip(z)*flip(z)),
   |z| <= 100
}

FlipProbJ1 { ; Ron Barnett
   ; try p1 = (1,1)
   z = pixel:
   z = flip(z)*(1-z) + p1,
   |z| <= 100
}

FlipProbJ2 { ; Ron Barnett
   ; try p1 = (-0.88,0.625)
   z = pixel:
   z = z*(p1-flip(z)) + p1,
   |z| <= 100
}

FlipProbM1 { ; Ron Barnett
   ; provides a "map" of locations for FlipProbJ1
   z = pixel:
   z = flip(z)*(1-z) + pixel,
   |z| <= 100
}

FlipProbM2 { ; Ron Barnett
   ; provides a "map" of locations for FlipProbJ2
   z = pixel:
   z = z*(pixel-flip(z)) + pixel,
   |z| <= 100
}

Fly(XAXIS_NOPARM) { ; Jonathan Osuch
   ; Generalized by Tobey J. E. Reed
   ; Try p1=0, p2=4, fn1=sqr, fn2=sqr
   ; Note:  use floating point
   z = p1:
   x = real(z),
   (x<0)*(z = fn1(z)+pixel),
   (0<=x)*(z = fn2(z)-pixel),
   |z|< = p2
}

FlyC(XAXIS_NOPARM) { ; Jonathan Osuch
   ; Generalized by Tobey J. E. Reed
   ; Try p1=0, p2=4, fn1=sqr, fn2=sqr
   ; Note:  use floating point
   z=p1:
   x=real(z),
   (z=fn1(z)+pixel)*(x<0)+(z=fn2(z)-pixel)*(0<=x),
   |z|<=p2
}

FlyingSquirrel(XAXIS_NOPARM) { ; Jonathan Osuch
   ; Generalized by Tobey J. E. Reed
   ; Try p1=0, p2=4, fn1=sin, fn2=cosxx, fn3=sqr
   ; Note:  use floating point
   z  = p1, x  = |z|:
   (1<x) * (z=fn1(z) / fn2(z) + pixel),
   z  = fn3(z)+pixel,
   x  = |z|,
   x <= p2
}

FlyingSquirrelC(XAXIS_NOPARM) { ; Jonathan Osuch
   ; Generalized by Tobey J. E. Reed
   ; Try p1=0, p2=4, fn1=sin, fn2=cos, fn3=sqr
   ; Note:  use floating point
   z=p1, x=|z|:
   (z=fn1(z)/fn2(z)+pixel)*(1<x)+(z=z)*(x<=1),
   z=fn3(z)+pixel, x=|z|,
   x<=p2
}

Form3 (XAXIS) { ; Peter Lewman's formulas for Fractint.
   z = Pixel, c = Pixel:
   z = c * z * ( p1 - z ),
   |z| < 4
}

Form4 (XAXIS) { ; Peter Lewman's formulas for Fractint.
   z = Pixel, c = P1:
   z = c * z * ( p2 - z ),
   |z| < 4
}

Form5 (XAXIS) { ; Peter Lewman's formulas for Fractint.
   z = Pixel, c = Pixel:
   z = p1 / ( fn1(z) + c ),
   |z| < 4
}

Form6 (XAXIS) { ; Peter Lewman's formulas for Fractint.
   z = Pixel, c = Pixel:
   z = z^6 + fn1(z) + c,
   |z| < 4
}

Form7 (XYAXIS) { ; Peter Lewman's formulas for Fractint.
   z = Pixel, c = Pixel:
   z = ( c * fn1( fn2(z) + 1 ) ) / ( z * ( fn3(z) - 1) ),
   |z| < 4
}

FractalFender(XAXIS_NOPARM) { ; Jonathan Osuch
   ; Generalized by Tobey J. E. Reed
   ; Try p1=0, p2=4, fn1=cosh, fn2=sqr
   ; Try p1=0, p2=4, fn1=cosxx, fn2=sqr
   ; Note:  use floating point
   z  = p1, x  = |z|:
   (1<x) * (z=fn1(z)+pixel),
   z  = fn2(z)+pixel,
   x  = |z|,
   x <= p2
}

Frame-RbtJ { ; Ron Barnett
   ; try p1 = (-1.37, 0.57)
   z = pixel:
   z = z*z*z/5 + z*z + p1,
   |z| <= 100
}

Frog(XAXIS_NOPARM) { ; Jonathan Osuch
   ; Generalized by Tobey J. E. Reed
   ; Try p1=0, p2=4, fn1=tanh, fn2=sqr
   ; Note:  use floating point
   z  = p1, x  = |z|:
   (1<x) * (z=fn1(z) + pixel),
   z  = fn2(z)+pixel,
   x  = |z|,
   x <= p2
}

FrogC(XAXIS_NOPARM) { ; Jonathan Osuch
   ; Generalized by Tobey J. E. Reed
   ; Try p1=0, p2=4, fn1=tanh, fn2=sqr
   ; Note:  use floating point
   z  = p1, x  = |z|:
   (z  = fn1(z)+pixel) * (1<x) + (z=z) * (x<=1),
   z  = fn2(z)+pixel,
   x  = |z|,
   x <= p2
}

FrRbtGenJ { ; Ron Barnett
   z = pixel:
   z = p1*z*z*z + z*z + p2,
   |z| <= 100
}

Fzpcocoh  { ; Lee Skinner
   z = pixel, f = 1. / cosh(pixel):
   z = cosxx (z)  + f,
   |z| <= 50
}

Fzpcopch  { ; Lee Skinner
   z = pixel, f = pixel ^ (cosh(pixel) ):
   z = cosxx (z)  + f,
   |z| <= 50
}

Fzpcopcs  { ; Lee Skinner
   z = pixel, f = pixel ^ (1. / cosxx(pixel) ):
   z = cosxx (z)  + f,
   |z| <= 50
}

Fzpcopct  { ; Lee Skinner
   z = pixel, f = pixel ^ (cosxx(pixel) / sin(pixel) ):
   z = cosxx (z)  + f,
   |z|<= 50
}

Fzpcophc  { ; Lee Skinner
   z = pixel, f = pixel ^ (1. / cosh(pixel) ):
   z = cosxx (z)  + f,
   |z| <= 50
}

Fzpcophs  { ; Lee Skinner
   z = pixel, f = pixel ^ (1. / sinh(pixel) ):
   z = cosxx (z)  + f,
   |z| <= 50
}

Fzpcopht  { ; Lee Skinner
   z = pixel, f = pixel ^ cotanh(pixel) :
   z = cosxx (z)  + f,
   |z|<= 50
}

Fzpcopse  { ; Lee Skinner
   z = pixel, f = pixel ^ (1. / sin(pixel) ):
   z = cosxx (z)  + f,
   |z| <= 50
}

Fzpcopsh  { ; Lee Skinner
   z = pixel, f = pixel ^ (sinh(pixel) ):
   z = cosxx (z)  + f,
   |z| <= 50
}

Fzpcopsq  { ; Lee Skinner
   z = pixel, f = pixel ^ (sqr(pixel) ):
   z = cosxx (z)  + f,
   |z| <= 50
}

Fzpcopta  { ; Lee Skinner
   z = pixel, f = pixel ^ (sin(pixel) / cosxx(pixel) ):
   z = cosxx (z)  + f,
   |z|<= 50
}

Fzpcopth  { ; Lee Skinner
   z = pixel, f = pixel ^ tanh(pixel):
   z = cosxx (z)  + f,
   |z|<= 50
}

Fzpcoseh  { ; Lee Skinner
   z = pixel, f = 1. / sinh(pixel):
   z = cosxx (z)  + f,
   |z| <= 50
}

Fzppchex  { ; Lee Skinner
   z = pixel, f = exp (pixel):
   z = cosh (z) + f,
   |z| <= 50
}

Fzppchlo  { ; Lee Skinner
   z = pixel, f = log (pixel):
   z = cosh (z) + f,
   |z| <= 50
}

Fzppchsh  { ; Lee Skinner
   z = pixel, f = sinh (pixel):
   z = cosh (z) + f,
   |z| <= 50
}

Fzppchsi  { ; Lee Skinner
   z = pixel, f = sin (pixel):
   z = cosh (z) + f,
   |z| <= 50
}

Fzppchsq  { ; Lee Skinner
   z = pixel, f = sqr (pixel):
   z = cosh (z) + f,
   |z| <= 50
}

Fzppcoch  { ; Lee Skinner
   z = pixel, f = cosh (pixel):
   z = cosxx (z)  + f,
   |z| <= 50
}

Fzppcocs  { ; Lee Skinner
   z = pixel, f = 1. / cosxx(pixel):
   z = cosxx (z)  + f,
   |z| <= 50
}

Fzppcoct  { ; Lee Skinner
   z = pixel, f = cosxx(pixel) / sin(pixel):
   z = cosxx (z)  + f,
   |z|<= 50
}

Fzppcoex  { ; Lee Skinner
   z = pixel, f = exp (pixel):
   z = cosxx (z)  + f,
   |z| <= 50
}

Fzppcohs  { ; Lee Skinner
   z = pixel, f = sinh (pixel):
   z = cosxx (z)  + f,
   |z| <= 50
}

Fzppcoht  { ; Lee Skinner
   z = pixel, f = cotanh(pixel):
   z = cosxx (z)+f,
   |z|<= 50
}

Fzppcolo  { ; Lee Skinner
   z = pixel, f = log (pixel):
   z = cosxx (z)  + f,
   |z| <= 50
}

Fzppcopc  { ; Lee Skinner
   z = pixel, f = pixel ^ (cosxx(pixel) ):
   z = cosxx (z)  + f,
   |z| <= 50
}

Fzppcope  { ; Lee Skinner
   z = pixel, f = pixel ^ (exp(pixel) ):
   z = cosxx (z)  + f,
   |z| <= 50
}

Fzppcopl  { ; Lee Skinner
   z = pixel, f = pixel ^ (log(pixel) ):
   z = cosxx (z)  + f,
   |z| <= 50
}

Fzppcopo  { ; Lee Skinner
   z = pixel, f = (pixel) ^ (pixel):
   z = cosxx (z)  + f,
   |z| <= 50
}

Fzppcopr  { ; Lee Skinner
   z = pixel, f = pixel ^ (1. / pixel):
   z = cosxx (z)  + f,
   |z| <= 50
}

Fzppcops  { ; Lee Skinner
   z = pixel, f = pixel ^ (sin(pixel) ):
   z = cosxx (z)  + f,
   |z| <= 50
}

Fzppcore  { ; Lee Skinner
   z = pixel, f = 1. / (pixel):
   z = cosxx (z)  + f,
   |z| <= 50
}

Fzppcos   { ; Lee Skinner
   z = pixel, f = cosxx (pixel):
   z = cosxx (z)  + f,
   |z| <= 50
}

Fzppcose  { ; Lee Skinner
   z = pixel, f = 1. / sin(pixel):
   z = cosxx (z)  + f,
   |z| <= 50
}

Fzppcosh  { ; Lee Skinner
   z = pixel, f = cosh (pixel):
   z = cosh (z) + f,
   |z| <= 50
}

Fzppcosi  { ; Lee Skinner
   z = pixel, f = sin (pixel):
   z = cosxx (z)  + f,
   |z| <= 50
}

Fzppcota  { ; Lee Skinner
   z = pixel, f = sin(pixel) / cosxx(pixel):
   z = cosxx (z)  + f,
   |z|<= 50
}

Fzppcoth  { ; Lee Skinner
   z = pixel, f = tanh(pixel):
   z = cosxx (z)+f,
   |z|<= 50
}

Fzppexch  { ; Lee Skinner
   z = pixel, f = cosh (pixel):
   z = exp (z)  + f,
   |z| <= 50
}

Fzppexco  { ; Lee Skinner
   z = pixel, f = cosxx (pixel):
   z = exp (z)  + f,
   |z| <= 50
}

Fzppexlo  { ; Lee Skinner
   z = pixel, f = log (pixel):
   z = exp (z)  + f,
   |z| <= 50
}

Fzppexp   { ; Lee Skinner
   z = pixel, f = exp (pixel):
   z = exp (z)  + f,
   |z| <= 50
}

Fzppexsh  { ; Lee Skinner
   z = pixel, f = sinh (pixel):
   z = exp (z)  + f,
   |z| <= 50
}

Fzppexsi  { ; Lee Skinner
   z = pixel, f = sin (pixel):
   z = exp (z)  + f,
   |z| <= 50
}

Fzppexsq  { ; Lee Skinner
   z = pixel, f = sqr (pixel):
   z = exp (z)  + f,
   |z| <= 50
}

Fzppshch  { ; Lee Skinner
   z = pixel, f = cosh (pixel):
   z = sinh (z) + f,
   |z| <= 50
}

Fzppshco  { ; Lee Skinner
   z = pixel, f = cosxx (pixel):
   z = sinh (z) + f,
   |z| <= 50
}

Fzppshex  { ; Lee Skinner
   z = pixel, f = exp (pixel):
   z = sinh (z) + f,
   |z| <= 50
}

Fzppshlo  { ; Lee Skinner
   z = pixel, f = log (pixel):
   z = sinh (z) + f,
   |z| <= 50
}

Fzppshsi  { ; Lee Skinner
   z = pixel, f = sin (pixel):
   z = sinh (z) + f,
   |z| <= 50
}

Fzppshsq  { ; Lee Skinner
   z = pixel, f = sqr (pixel):
   z = sinh (z) + f,
   |z| <= 50
}

Fzppsich  { ; Lee Skinner
   z = pixel, f = cosh (pixel):
   z = sin (z)  + f,
   |z| <= 50
}

Fzppsico  { ; Lee Skinner
   z = pixel, f = cosxx (pixel):
   z = sin (z)  + f,
   |z| <= 50
}

Fzppsiex  { ; Lee Skinner
   z = pixel, f = exp (pixel):
   z = sin (z)  + f,
   |z| <= 50
}

Fzppsinh  { ; Lee Skinner
   z = pixel, f = sinh (pixel):
   z = sinh (z) + f,
   |z| <= 50
}

Fzppsish  { ; Lee Skinner
   z = pixel, f = sinh (pixel):
   z = sin (z)  + f,
   |z| <= 50
}

Fzppsisq  { ; Lee Skinner
   z = pixel, f = sqr (pixel):
   z = sin (z)  + f,
   |z| <= 50
}

Fzppsqlo  { ; Lee Skinner
   z = pixel, f = log (pixel):
   z = sqr (z)  + f,
   |z| <= 50
}

Fzppsqsh  { ; Lee Skinner
   z = pixel, f = sinh (pixel):
   z = sqr (z)  + f,
   |z| <= 50
}

Fzppsqsi  { ; Lee Skinner
   z = pixel, f = sin (pixel):
   z = sqr (z)  + f,
   |z| <= 50
}

GLYNN(XAXIS) { ; Based on an illustration in Science PROBE!  and a
   ; formula by Earl Glynn in Computers and the Imagination,
   ; by Clifford Pickover.   Try p1 = 1.5, p2 = -0.2
   ; Jon Horner, FRAC'Cetera
   ;
   z = pixel :
   z = z ^ p1 + p2 ,
   |z| <=4
}

Gopalsamy1 { ; Ron Barnett
   ; try p1 = (0.29,0.29)
   z = pixel:
   x = real(z), y = imag(z),
   x1 = -2*x*y + p1,
   y = y*y - x*x,
   z = x1 + flip(y),
   |z| <= 4
}

Gopalsamy2 { ; Ron Barnett
   ; try p1 = 0.25
   z = pixel:
   x = real(z), y = imag(z),
   x1 = -4*x*y + p1,
   y = 4*y*y - x*x,
   z = x1 + flip(y),
   |z| <= 4
}

Gopalsamy3 { ; Ron Barnett
   ; try p1 = 1.099
   z = pixel:
   x = real(z), y = imag(z),
   x1 = 3*x*y*y - x*x*x + p1,
   y = y*y*y - 3*x*x*y,
   z = x1 + flip(y),
   |z| <= 4
}

Gopalsamy4 { ; Ron Barnett
   ; p1 = 0.31
   z = pixel:
   x = real(z), y = imag(z),
   x1 = -x*y + p1,
   y = 2*y*y - 3*x*x,
   z = x1 + flip(y),
   |z| <= 4
}

Gopalsamy5 { ; Ron Barnett
   ; try p1 = 0.835
   z = pixel:
   x = real(z), y = imag(z),
   x1 = 2*x*y,
   y1 = x*x - y*y,
   x2 = -2*x1*y1 + p1,
   y = y1*y1 - x1*x1,
   z = x2 + flip(y),
   |z| <= 4
}

GopalsamySin2 { ; Ron Barnett
   ; use floating point
   z = pixel:
   x = real(z), y = imag(z),
   x1 = sin(x)*cosh(y),
   y1 = cos(x)*sinh(y),
   x2 = -2*x1*y1 + p1,
   y = y1*y1 - x1*x1,
   z = x2 + flip(y),
   |z| <= 100
}

GopalsamySin { ; Ron Barnett
   z = pixel:
   x = real(z), y = imag(z),
   x1 = -sin(x)*cosh(y) + p1,
   y = -cos(x)*sinh(y),
   z = x1 + flip(y),
   |z| <= 100
}

GopalsamyExp { ; Ron Barnett
   z = pixel:
   x = real(z), y = imag(z),
   x1 = -exp(x)*cos(y) + p1,
   y = -exp(x)*sin(y),
   z = x1 + flip(y),
   |z| <= 100
}

GopalsamyExp2 { ; Ron Barnett
   z = pixel:
   x = real(z), y = imag(z),
   x1 = exp(x)*cos(y),
   y1 = exp(x)*sin(y),
   x2 = -2*x1*y1 + p1,
   y = y1*y1 - x1*x1,
   z = x2 + flip(y),
   |z| <= 100
}

GopalsamySinh2 { ; Ron Barnett
   z = pixel:
   x = real(z), y = imag(z),
   x1 = sinh(x)*cos(y),
   y1 = cosh(x)*sin(y),
   x2 = -2*x1*y1 + p1,
   y = y1*y1 - x1*x1,
   z = x2 + flip(y),
   |z| <= 100
}

GopalsamySinh { ; Ron Barnett
   z = pixel:
   x = real(z), y = imag(z),
   x1 = -sin(x)*cosh(y) + p1,
   y = -cos(x)*sinh(y),
   z = x1 + flip(y),
   |z| <= 100
}

IfThenfn1fn2(XAXIS_NOPARM) { ; Jonathan Osuch
   ; Generalized by Tobey J. E. Reed
   ; Try p1=0, p2=4, fn1=sin, fn2=cos
   ; Note:  use floating point
   z  = p1, x  = |z|:
   (z  = fn1(z)) * (1<x)+(z=z)*(x<=1),
   (z  = fn2(z)  + pixel),
   x  = |z|,
   x <= p2
}

IfThenElsefn1fn2(XAXIS_NOPARM) { ; Jonathan Osuch
   ; Generalized by Tobey J. E. Reed
   ; Try p1=0, p2=4, fn1=sqr, fn2=sin
   ; Note:  use floating point
   z  = p1, x  = |z|:
   (z  = fn1(z)+pixel) * (1<x)+(z=fn2(z)+pixel) * (x<=1),
   x  = |z|,
   x <= p2
}

IfElsefn1fn2fn3(XAXIS_NOPARM) { ; Jonathan Osuch
   ; Generalized by Tobey J. E. Reed
   ; Try p1=0, p2=4, fn1,2,3=whatever
   ; Note:  use floating point
   z   = p1, x   = 1:
   (z=fn1(z)+pixel)*(x<10)+(z=fn2(z)+pixel)*(10<=x)*(x<20)+(z=fn3(z)+pixel)*(20<=x),
   x   = x+1,
   |z| <= p2
}

IkeFrRbtGenJ { ; Ron Barnett
   z = pixel:
   z = p1*z*z*z + (p2-1)*z*z - p2,
   |z| <= 100
}

IkeFrRbtGenM { ; Ron Barnett
   z = 2*(1-pixel)/(3*p1):
   z = p1*z*z*z + (pixel-1)*z*z - pixel,
   |z| <= 100
}

IkeGenJ { ; Ron Barnett
   z = pixel:
   z =p1*z*z*z + (p2-1)*z - p2,
   |z| <= 100
}

IkeGenM { ; Ron Barnett
   z = ((1-pixel)/(3*p1))^0.5:
   z =p1*z*z*z + (pixel-1)*z - pixel,
   |z| <= 100
}
