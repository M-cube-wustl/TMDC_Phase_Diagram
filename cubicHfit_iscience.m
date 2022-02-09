close all, fclose all
clear

dir = 'H_data';

isQuart = false;
isXtraTab = false;
isChalc = true;
fps = fps_in_dir(dir);
imdir = dircat(dir,'cubicfit');
dldir = dircat(dir,'datalines');
figsize = [5,5];
mkdir(imdir);
mkdir(dldir);
N = numel(fps);
%H = zeros(5,N);
names = cell(1,N);
dx = .01;
if isXtraTab
    shift = -1;
else
    shift = 0;
end
for i=1:N
    fp = fps{i};
    N_col = columnsInFile(fp)+shift;
    fid = fopen(fp,'r');
    fn = file_from_fp(fps{i});
    input = fscanf(fid,'%f',[N_col 2])';
    fclose(fid);
    M = size(input,2);
    Hin(1:M,i) = 1000*input(2,:)';
    X(1:M,i) = input(1,:)';
    names{i} = fn;
end

k = 8.61733*10^-5;
if isChalc
    k = 2*k;
end    
close all

namestr = strings(1,N);
headers = ["a","b","max","asym","Tmisc"];
D = zeros(N,numel(headers));

for i=1:N
    
    
    name = names{i};
    alloy = name;
    phase = '2H';
    figure;
    hold on
    set(gca,'fontsize', 18);
    imgfp = dircat(imdir,[name,'-H.png']);
    
    xlabel('$x$','Interpreter','latex')
    ylabel('$\Delta H $ (meV)','Interpreter','latex')
    y = Hin(:,i);
    H = y';
    
    xfun = @(x) [x.^3-x,x-x.^2];
    xfun4 = @(x) [x.^4-x,x.^3-x,x.^2-x];
    x = X(:,i);
    xs = (0:dx:1)';

    d1 = @(x,B) B(1)*(3*x.^2-1)+B(2)*(1-2*x);
    d2 = @(x,B) B(1)*6*x-B(2)*2;
    if isQuart
        d2 = @(x,B) 2*(3*x.*(2*B(1)*x+B(2))+B(3));
    end
    
    [ B3,r2_3 ] = multilinearreg( x, y, xfun );
    [ B4,r2_4 ] = multilinearreg( x, y, xfun4 );
    
    H1 = @(x) xfun(x)*B3;
    H2 = @(x) xfun(x)*B3;
    
    if isQuart
        H1 = @(x) xfun4(x)*B4;
        H2 = @(x) xfun4(x)*B4;
    end
    nperfu = 1;
    if isChalc
        nperfu = 2;
    end
    a = B3(1);
    b = 3/2*B3(1) - B3(2);
    
    D(i,1) = a;
    D(i,2) = b;
    if b<0
        D(i,3) = ((3*a-2*b)-sqrt((2*b-3*a)^2-12*a*(a/2-b)))/(6*a)-.5;
    else
        D(i,3) = ((3*a-2*b)+sqrt((2*b-3*a)^2-12*a*(a/2-b)))/(6*a)-.5;
    end
    D(i,4) = sign(a)*sign(b)*abs(a)^(1/3)/abs(b)^(1/2);
    xspin = 0:.001:1;
    spin = @(x) -1/k*(3*a*(2*x-1)+2*b).*x.*(1-x)/1000;
    if isQuart
        spin = @(x) -1/k*2*(3*x.*(2*B4(1)*x+B4(2))+B4(3)).*x.*(1-x)/1000;
    end
    [Tmax,imax] = max(spin(xspin));
    [Tmin, imin] = min(spin(xspin));
    if abs(Tmax)>abs(Tmin)
        Tmisc = Tmax;
        x0 = xspin(imax);
    else
        Tmisc = Tmin;
        x0 = xspin(imin);
    end
    if mean(H)<0
        Tmisc = 0
    end
        
    D(i,5) = Tmisc;
    
    H_fit = H1(xs);
    xx = xs;
    %plot(xs,H1(xs),'LineWidth',2.5)
    %plot(x,y,'o','markers',15)
    ylim([-20,100])
    ylims = ylim;
    xlim([0 1])
    if y(3)>0
        plot_vert(B3(2)/(3*B3(1)),ylims(1),ylims(2),2.5)
    end
    
    r = r2_4;
    t = strcat(name,' - $r^2=$',num2str(r2_4));
    title(t,'Interpreter','latex')
    set(gcf, 'PaperSize', figsize);
    %saveas(gcf,imgfp);
    

    
    %s = spline(x,y,xs);
    %lap = del2(s,dx);
    
    
    %plot(xs,lap)
    
    spin = @(x,ddH) ddH/k.*x.*(x-1)/1000;
    %mspin = @(x) -spin(x,d2(xs,B3));
    %opts = optimset('Display','iter');
    %Tmisc = fminsearch(mspin,.5);
    %D(i,5) = Tmisc;
    xaxis = zeros(numel(xs),1);
    figure
    imgfp = dircat(imdir,[name,'-phase.png']);
    set(gca,'fontsize', 18);
    set(gcf, 'PaperSize', figsize);
    hold on
    xS=xs;
    if isQuart
        yS=spin(xs,d2(xs,B4))';
        fill([xs',xs'],[xaxis',yS],'b','LineWidth',2.5)
        spinfun = @(x) spin(x,d2(x,B4));
    else
        yS=spin(xs,d2(xs,B3))';
        fill([xs',xs'],[xaxis',yS],'b','LineWidth',2.5)
        spinfun = @(x) spin(x,d2(x,B3));
    end

    %plot(xs,spin(xs,lap))
    xlabel('$x$','Interpreter','latex')
    ylabel('Temperature (K)','Interpreter','latex')
    ylims = [0,1200];
    %ylim([0 ylims(2)])
    t = strcat('Phase Diagram - ',name);
    title(t,'Interpreter','latex')
    ylim(ylims)
    
    saveas(gcf,imgfp);
    
    Tmisc
    Tset = [150:25:(Tmisc*.7)];
    Tset = [150:10:(Tmisc*.99)];
    X0 = [.001,.999];
    %if strcmp(name,'Mo$_{1-x}$V$_x$S$_2$')
    
    Tiso=0;
    
        %else
        %Tiso = 0;
    %end
    
    
    if isQuart
        spinfun = @(x) spin(x,d2(x,B4))-Tiso;
    else
        spinfun = @(x) spin(x,d2(x,B3))-Tiso;
    end
    %options = optimoptions('lsqnonlin','StepTolerance',StepTolerance,'MaxFunctionEvaluations',MaxFunctionEvaluations,'FunctionTolerance',FunctionTolerance);
    if mean(H)>0
        [x_spin1] = lsqnonlin(spinfun,.1,[0 0],[1,1]);%,options);
        [x_spin2] = lsqnonlin(spinfun,.9,[0 0],[1,1]);%,options);
    else
        x_spin1=0;
        x_spin2=0;
    end
    if numel(Tset)>0
        x_spin = [x_spin1, x_spin2];
    else
        x_spin = 0;
    end
    
    [ xB,X2,Tset,x_bin,y_bin,xG,yG,y_spin] = commontangent_fun( H1,H2,X0,nperfu,Tset,x0,Tmisc,Tiso,x_spin );
    
    labels = {'alloy', 'phase', 'xS', 'yS', 'xB', 'yB', 'H', 'x', 'xx', 'H_fit','x_spin','y_spin','x_bin','y_bin','xG','yG','Tiso'};
    datae = {alloy, phase, xS, yS, xB, Tset, H/1000, x, xx, H_fit/1000,x_spin,y_spin,x_bin,y_bin,xG,yG,Tiso};
    dlfp = dircat(dldir,[name,'.dl']);

    write_DataLines( dlfp, labels, datae )

end
write_ALL_DATA( imdir, 'fitdata', names, headers, D, 'fitdata' )