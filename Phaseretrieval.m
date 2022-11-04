clear all;clc;close all;
lambda=1.2*10^-9;
k=2*pi/lambda;

z=20;
d=0.2;

dx2=4*10^-6;
M=512;
L2=M*dx2;

dx1=lambda*z/M/dx2; %image sample interval
L1=lambda*z/dx2; %obs sidelength

m = 128; % size of the sample image��Ʒ��С
m2 = m/2;

sample = rgb2gray(imread('Fudan.jpg')); 
sample = imresize(sample,[m m]);
sample = MatMap(sample,0,1); % the sample matrix is normalized

figure;imshow(sample,'InitialMagnification',200);axis on;title('Intensity of the object domain');set(gca,'xticklabel',{'a','b','c'});

S = zeros(M,M);
S(M/2-(m2-1):M/2+1+(m2-1) ,M/2-(m2-1):M/2+1+(m2-1)) = sample; %��sample����S���м�

 sup = circle_mask(M,m,M/2,M/2); % generate a circle support
%  sup = triMask(M,m/2+8,M/2+10,M/2); % Ҫ�ǻ���һ����������ģЧ�������

S = S.*sup;%����Ʒ����Ĥ
figure;imshow(S);axis on;title('֧���������µ��﷽ǿ�ȷֲ�');

S = abs(fftshift(fft2(S))); % generate the modulus of the diffraction pattern����Ҷ�任��������ͼ���ķ�ֵ
%figure;imshow(S);axis square;title('The modulus of diffraction pattern');
figure;
imagesc(nthroot(S.*S,3));
axis square;
title('The diffraction pattern');
%%

itnum = 500; % iteration number
g = rand(M,M);%g���﷽�������������񷽳����渵��Ҷ�任

figure;
for i = 1:itnum
 
    %=================ER========================
 if mod(floor(i/50),8)==0
        g = projectSup(projectM(g,S),sup);
 end
    %=================HIO========================
if mod(floor(i/50),8)~=0
     g2 = projectM(g,S);
     g3 = g2.*sup;
     g = (g3>0).*g2 + (g3<=0).*(g-0.99.*g2);
end
    %==============display the reconstruct sample image===================
  %  subplot(2,1,1);
    imshow(g(M/2-(m2-1):M/2+1+(m2-1),M/2-(m2-1):M/2+1+(m2-1)),'InitialMagnification',200);
    title(strcat('�ؽ��������������',num2str(i)));
    pause(0.01); % ÿһ�������������ʾʱ��
    hold on;
    % ��ʾ�м�ĵ������ʱֻ��ʾ N x N ��С�ؽ�ͼ���м�� m x m ��С�������Աߴ�Ƭ�ĺ�ɫ���ֲ���ʾ
    
    Ef(i,1)=(zeros(1,M)+1)*(abs(fftshift(fft2(g)))-S).^2*(zeros(M,1)+1)/M^2;
    Eo(i,1)=(zeros(1,M)+1)*(g.*(1-sup)).^2*(zeros(M,1)+1);
   
end
 figure;plot(Ef);
 figure;plot(Eo);

    

