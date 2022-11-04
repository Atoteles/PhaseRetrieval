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

m = 128; % size of the sample image样品大小
m2 = m/2;

sample = rgb2gray(imread('Fudan.jpg')); 
sample = imresize(sample,[m m]);
sample = MatMap(sample,0,1); % the sample matrix is normalized

figure;imshow(sample,'InitialMagnification',200);axis on;title('Intensity of the object domain');set(gca,'xticklabel',{'a','b','c'});

S = zeros(M,M);
S(M/2-(m2-1):M/2+1+(m2-1) ,M/2-(m2-1):M/2+1+(m2-1)) = sample; %把sample放在S的中间

 sup = circle_mask(M,m,M/2,M/2); % generate a circle support
%  sup = triMask(M,m/2+8,M/2+10,M/2); % 要是换成一个三角形掩模效果会更好

S = S.*sup;%给样品加掩膜
figure;imshow(S);axis on;title('支持域作用下的物方强度分布');

S = abs(fftshift(fft2(S))); % generate the modulus of the diffraction pattern傅里叶变换生成衍射图样的幅值
%figure;imshow(S);axis square;title('The modulus of diffraction pattern');
figure;
imagesc(nthroot(S.*S,3));
axis square;
title('The diffraction pattern');
%%

itnum = 500; % iteration number
g = rand(M,M);%g是物方场，理论上是像方场的逆傅里叶变换

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
    title(strcat('重建结果：迭代步数',num2str(i)));
    pause(0.01); % 每一步迭代结果的显示时间
    hold on;
    % 显示中间的迭代结果时只显示 N x N 大小重建图像中间的 m x m 大小的区域，旁边大片的黑色部分不显示
    
    Ef(i,1)=(zeros(1,M)+1)*(abs(fftshift(fft2(g)))-S).^2*(zeros(M,1)+1)/M^2;
    Eo(i,1)=(zeros(1,M)+1)*(g.*(1-sup)).^2*(zeros(M,1)+1);
   
end
 figure;plot(Ef);
 figure;plot(Eo);

    

