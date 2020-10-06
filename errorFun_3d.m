function fun=errorFun_3d(x)
% global Pot_Mat len;
% %scatter3(x(:,1),x(:,2),x(:,3),5,1:3:3*(len+1),'filled','MarkerFaceAlpha',0.01)
% alpha = x(len+1);
% A_mat = repmat(x(1:len,4),1,len); %replaced x(1:len,4) with A_vect
% v0 = x(len+1,2);
% r = dist(x(1:len,1:3)');
% prefun1 = (log(Pot_Mat-v0) - log(A_mat) + alpha*log(r));
% prefun1(eye( size( prefun1 ) )==1) = 0; %zeros all the terms in the diagonal which were zero before and became inf after log
% prefun2 = sum(sum(   (prefun1).^6   )) ;
% fun = prefun2;
global Pot_Mat len;
alpha = x(len+1);
v0 = x(len+1,2);
r = dist(x(1:len,1:3)');
prefun = zeros(len);
for idx = 1:len
    for jdx = 1:len
        if idx == jdx
            continue
        end
        prefun(idx,jdx) = log(abs(Pot_Mat(idx,jdx))-v0) - log(x(jdx,4)) + alpha*log(r(idx,jdx));
    end
end
fun = sum(sum(prefun.^2));
end
