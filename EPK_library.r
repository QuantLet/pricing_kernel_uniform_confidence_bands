###########################################################

K = function (u)
{
    kk=length(u);
    ret=vector(, length=kk);
    for (i in c(1:kk))
    {
        ret[i]=dnorm(u[i]);         # Gaussian
#       if (abs(u[i])>1) {ret[i]=0;} else {ret[i]=(3*(1-u[i]^2)/4) };       # Epanechnikov
    }
    return(ret);
}

K.h = function (u, h)
{
    kk=length(u);
    ret=vector(, length=kk);
    for (i in c(1:kk)){ret[i]=K(u[i]/h)/h;}
    return(ret);
}

###########################################################

rookley = function(money, sigma, sigma1, sigma2, r, tau) 
{
  rm=length(money);
  sqrttau=sqrt(tau);
  exprt=exp(r*tau);
  rtau=r*tau;

    d1=(log(money)+tau*(r+0.5*sigma^2))/(sigma*sqrttau);
    d2=d1-sigma*sqrttau;
    
    d11 = 1/(money*sigma*sqrttau)+(-(log(money)+tau*r)/(sqrttau*sigma^2)+0.5*sqrttau)*sigma1;
        
    d21 = 1/(money*sigma*sqrttau)+(-(log(money)+tau*r)/(sqrttau*sigma^2)-0.5*sqrttau)*sigma1;

    d12 = -1/(sqrttau*money^2*sigma)-sigma1/(sqrttau*money*sigma^2)+ sigma2*(0.5*sqrttau-(log(money)+r*tau)/(sqrttau*sigma^2))+ sigma1*(2*sigma1*(log(money)+r*tau)/(sqrttau*sigma^3)-1/(sqrttau*money*sigma^2));
    d22 = -1/(sqrttau*(money^2)*sigma)-sigma1/(sqrttau*money*sigma^2) + sigma2*(-0.5*sqrttau-(log(money)+r*tau)/(sqrttau*sigma^2))+ sigma1*(2*sigma1*(log(money)+r*tau)/(sqrttau*sigma^3)-1/(sqrttau*money*sigma^2));

    f = pnorm(d1)-pnorm(d2)/(money*exprt);
    f1=dnorm(d1)*d11-dnorm(d2)*d21/(exprt*money)+pnorm(d2)/(exprt*money^2);
    f2=dnorm(d1)*d12-d1*dnorm(d1)*d11^2-dnorm(d2)*d22/(exprt*money)+dnorm(d2)*d21/(exprt*money^2) + d2*dnorm(d2)*d21^2/(exprt*money)-2*pnorm(d2)/(exprt*money^3)+dnorm(d2)*d21/(exprt*money^2);

#cat("d1 ", d1[ngrid]," d2 ",d2[ngrid]," d11 ", d11[ngrid]," d12 ", d12[ngrid]," d21 ",d21[ngrid]," d22 ",d22[ngrid])
    return(cbind(f, f1, f2));
}

#########################       kernel functions
K = function (u)
{
    kk=length(u);
    ret=vector(, length=kk);
    for (i in c(1:kk))
    {
        ret[i]=dnorm(u[i]);         # Gaussian
#       if (abs(u[i])>1) {ret[i]=0;} else {ret[i]=(3*(1-u[i]^2)/4) };       # Epanechnikov
    }
    return(ret);
}

K.h = function (u, h)
{
    kk=length(u);
    ret=vector(, length=kk);
    for (i in c(1:kk)){ret[i]=K(u[i]/h)/h;}
    return(ret);
}


K.prime = function (u)
{
    kk=length(u);
    ret=vector(, length=kk);
    for (i in c(1:kk))
    {
        ret[i]=-u[i]*dnorm(u[i]);               # Gaussian
#       if(abs(u[i])>1) {ret[i]=0} else {ret[i]=(-3*u[i]/2)};       #   Epanechnikov
    }
    return(ret);
}

Kjp = function (u, j, N)
{
        kk=length(u);
        ret=vector(, length=kk);
        
        for (i in c(1:kk))
        {
            M = N;
            M[,j] = u[i]^{c(0:(dim(N)[1]-1))};
            if(abs(u[i])>1) {ret[i]=0;} else {ret[i]= K(u[i])*det(M)/det(N); };
        }
        return(ret);
}



######################################## N, Q, C matrices

#p=2;h=0.65; alpha=0.05;
func.L = function(p, h, alpha)
{
    N = matrix(0, nrow=p+1, ncol=p+1);
    Q = matrix(0, nrow=p+1, ncol=p+1);
    C = vector(, length=p+1);
    L = vector(, length=p+1);

    for(i in c(0:p))
    {
        for(j in c(0:p))
        {
            if ((i>1) && (j>1)) {Q[i+1,j+1]= integrate(function (xx) xx^(i+j)*(K.prime(xx))^2, -Inf,Inf)$value- 0.5 *(i*(i-1)+j*(j-1))*  integrate(function (x) x^(i+j-2)*(K(x))^2, -Inf,Inf)$value;}
            if ((i<=1) && (j<=1)) {Q[i+1,j+1]= integrate(function (xx) xx^(i+j)*(K.prime(xx))^2, -Inf,Inf)$value;}
            if ((i<=1) && (j>1)) {Q[i+1,j+1]= integrate(function (xx) xx^(i+j)*(K.prime(xx))^2, -Inf,Inf)$value- 0.5 *j*(j-1)*  integrate(function (x) x^(i+j-2)*(K(x))^2, -Inf,Inf)$value;}
            if ((i>1) && (j<=1)) {Q[i+1,j+1]= integrate(function (xx) xx^(i+j)*(K.prime(xx))^2, -Inf,Inf)$value- 0.5 *i*(i-1)*  integrate(function (x) x^(i+j-2)*(K(x))^2, -Inf,Inf)$value;}
            N[i+1,j+1]= integrate(function (xx) xx^(i+j)*K(xx), -1,1)$value;
        }
    }

    N.i = solve(N);

    for(j in c(0:p))
    { 
        C[j+1]= (N.i%*%Q%*%N.i)[j+1,j+1] /  integrate(function (xx) (Kjp(xx, j, N))^2, -Inf,Inf)$value;    
        L[j+1] = factorial(j)* (n*h^(2*j+1))^(-0.5) * ( (-2*log(h))^(1/2) + (-2*log(h))^(-1/2)*(-log(-0.5*log(1-alpha)) + log(sqrt(C[j+1])/(2*pi))));
    }
    print(C);
    return(L);
}

######################################### asymptotic bands


bands = function(p, h, hat.theta, x, y, vx)
{
    L = func.L(p, h, 0.05);sigma2=1;
    up.band= matrix(0,nrow=length(vx), ncol=p+1);
    lo.band= matrix(0,nrow=length(vx), ncol=p+1);
    V.all = matrix(0,nrow=length(vx), ncol=p+1);

    for(i  in c(1:length(vx)))
    {
        xx=vx[i];
        B.V = matrix(0, nrow=p+1, ncol=p+1);
        K.V = matrix(0, nrow=p+1, ncol=p+1);
    
        for(j in c(1:length(x)))
        {
#cat("ditheta", dim(hat.theta),"length(vx)",length(vx),"length(vx)",length(x),"length(y)", length(y), "length(x)", length(x),"i",i,"j",j,"\n");
            B.V = B.V + K.h(x[j]-xx,h) * (1/sigma2) * (solve(diag(h^c(0:p))) %*% t(t((x[j]-xx)^c(0:p))) %*% t((x[j]-xx)^c(0:p)) %*% solve(diag(h^c(0:p))));
            K.V = K.V + h* ( K.h(x[j]-xx,h) )^2 * ((y[j]-sum(hat.theta[i,]*((x[j]-xx)^c(0:p))))/sigma2)^2  * (solve(diag(h^c(0:p))) %*% t(t((x[j]-xx)^c(0:p))) %*% t((x[j]-xx)^c(0:p)) %*% solve(diag(h^c(0:p))));
        }
        B.V=B.V/n; K.V=K.V/n;
    
        B.i = solve(B.V);
        V = B.i%*%K.V%*% B.i;
        V.all[i,]=diag(V);
            
        for(ip in c(1:(p+1)))
        {
            lo.band[i,ip] =  - L[ip] * sqrt(V[ip, ip]);
            up.band[i,ip] =  L[ip] * sqrt(V[ip, ip]);
        }
    }   
    mylist=list(L=L, V=V.all, lo.band=lo.band, up.band=up.band);
    return(mylist);
}

