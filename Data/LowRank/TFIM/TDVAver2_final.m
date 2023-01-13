clear all
spin_num = 3;
layer = 2;
R = 8;
Jz = 1;
h = 0.25;
gamma = 1;
avalue = zeros([1,R],'single');
avalue(1) = 1.;
bvalue = zeros([1,(2*spin_num-1)*layer],'single');
bvalue(1:spin_num) = pi;
para_initial = [avalue, bvalue];
L = generate_Liou_TFIM(spin_num, Jz, h, gamma);

%[U, U_diff_b] = generate_hardware_efficient_ansatz(spin_num, layer, bvalue);
%[rho, rho_diff_a, rho_diff_b] = rho_state(U, U_diff_b, R, spin_num, avalue, bvalue);
% M = matrix_elements(rho_diff_a, rho_diff_b, R);
% V = vec_elements(rho, rho_diff_a, rho_diff_b, L, R);
% mz = Average_Magnetization(U,spin_num,R,avalue);

[tlist, mz, mx, paras] = VQD(L, R, spin_num, para_initial, 0.002, 5000, layer);
save(strcat('N',num2str(spin_num),'L',num2str(layer),'R',num2str(R),'Jz',num2str(Jz),'h',num2str(h),'.mat'));
figure(1)
plot(tlist,mz,'b.');
saveas(gcf,strcat('mz','N',num2str(spin_num),'L',num2str(layer),'R',num2str(R),'Jz',num2str(Jz),'h',num2str(h),'.jpg'));
figure(2)
plot(tlist,mx,'b.');
saveas(gcf,strcat('mx','N',num2str(spin_num),'L',num2str(layer),'R',num2str(R),'Jz',num2str(Jz),'h',num2str(h),'.jpg'));

function r = rx(a)
    r = [[cos(a/2), -1i*sin(a/2)];[-1i*sin(a/2), cos(a/2)]];
end

function r = ry(a)
    r = [[cos(a/2),-sin(a/2)];[sin(a/2),cos(a/2)]];
end

function r = rzz(a)
    r = complex(eye(4,'single'));
    r(1,1) = exp(-1j*a/2);
    r(2,2) = exp(1j*a/2);
    r(3,3) = exp(1j*a/2);
    r(4,4) = exp(-1j*a/2);
end

function x = sigmax_multi(spin_num,index)
    sigmax = [[0,1];[1,0]];
    x = kron(eye(2^(index-1),'single'),sigmax);
    x = kron(x,eye(2^(spin_num-index),'single'));

end

function y = sigmay_multi(spin_num,index)
    sigmay = [[0,-1j];[1j,0]];
    y = kron(eye(2^(index-1),'single'),sigmay);
    y = kron(y,eye(2^(spin_num-index),'single'));
end

function z = sigmaz_multi(spin_num,index)
    sigmaz = [[1,0];[0,-1]];
    z = kron(eye(2^(index-1),'single'),sigmaz);
    z = kron(z,eye(2^(spin_num-index),'single'));
end

function zz = sigmaz_sigmaz_multi(spin_num,index)
    sigmaz =  [[1,0];[0,-1]];
    if index == spin_num
        zz = kron(sigmaz,eye(2^(spin_num-2),'single'));
        zz = kron(zz,sigmaz);
    else
        zz = kron(eye(2^(index-1),'single'),kron(sigmaz,sigmaz));
        zz = kron(zz,eye(2^(spin_num-index-1),'single'));
    end
end

function r = rzz_multi(spin_num,index,a)
    r = kron(eye(2^(index-1),'single'),rzz(a));
    r = kron(r,eye(2^(spin_num-index-1),'single'));
end

function cz = controlz(spin_num,index)
    conz = eye(4,'single');
    conz(4,4) = -1;
    cz = kron(eye(2^(index-1),'single'),conz);
    cz = kron(cz,eye(2^(spin_num-index-1),'single'));
end

function [U, U_diff_b] = generate_hardware_efficient_ansatz(spin_num, layer, bvalue)
    Rxs = complex(zeros([spin_num*layer,2,2],'single'));
    Rys = complex(zeros([spin_num*layer,2,2],'single'));
    czs = complex(zeros([spin_num-1,2^spin_num,2^spin_num],'single'));
    for i = 1:layer
        for j = 1:spin_num
            Rxs(spin_num*(i-1)+j,:,:) = rx(bvalue(2*spin_num*(i-1)+j));
            Rys(spin_num*(i-1)+j,:,:) = ry(bvalue(2*spin_num*(i-1)+j+spin_num));
        end
    end
    for i = 1:spin_num-1
        czs(i,:,:) = controlz(spin_num,i);
    end
    Rx = complex(zeros([layer,2^spin_num,2^spin_num],'single'));
    Ry = complex(zeros([layer,2^spin_num,2^spin_num],'single'));
    if spin_num>2
        cz1 = reshape(czs(1,:,:),2^spin_num,2^spin_num);
        cz2 = reshape(czs(2,:,:),2^spin_num,2^spin_num);
        for i = 1:layer
            x = reshape(Rxs(spin_num*(i-1)+1,:,:),2,2);
            y = reshape(Rys(spin_num*(i-1)+1,:,:),2,2);
            for j = 2:spin_num
                x = kron(x,reshape(Rxs(spin_num*(i-1)+j,:,:),2,2));
                y = kron(y,reshape(Rys(spin_num*(i-1)+j,:,:),2,2));
            end
            Rx(i,:,:) = x;
            Ry(i,:,:) = y;
        end
        for j = 3:spin_num-1
            if mod(j,2) == 1
                cz1 = cz1*reshape(czs(j,:,:),2^spin_num,2^spin_num);
            else
                cz2 = cz2*reshape(czs(j,:,:),2^spin_num,2^spin_num);
            end
        end
    end
    if spin_num == 2
        cz1 = reshape(czs(1,:,:),2^spin_num,2^spin_num);
        cz2 = eye(4,'single');
        for i = 1:layer
            x = reshape(Rxs(spin_num*(i-1)+1,:,:),2,2);
            y = reshape(Rys(spin_num*(i-1)+1,:,:),2,2);
            for j = 2:spin_num
                x = kron(x,reshape(Rxs(spin_num*(i-1)+j,:,:),2,2));
                y = kron(y,reshape(Rys(spin_num*(i-1)+j,:,:),2,2));
            end
            Rx(i,:,:) = x;
            Ry(i,:,:) = y;
        end
    end
    cz = cz2*cz1;

    ULs = complex(zeros([2*layer,2^spin_num,2^spin_num],'single'));
    ULs(1,:,:) = reshape(Rx(1,:,:),2^spin_num,2^spin_num);
    ULs(2,:,:) = reshape(Ry(1,:,:),2^spin_num,2^spin_num)*reshape(ULs(1,:,:),2^spin_num,2^spin_num);
    ULs(2,:,:) = cz*reshape(ULs(2,:,:),2^spin_num,2^spin_num);
%     U = eye(2^spin_num,'single');
    for i = 2:layer
        ULs(2*i-1,:,:) = reshape(Rx(i,:,:),2^spin_num,2^spin_num)*reshape(ULs(2*i-2,:,:),2^spin_num,2^spin_num);
        ULs(2*i,:,:) = cz*reshape(Ry(i,:,:),2^spin_num,2^spin_num)*reshape(ULs(2*i-1,:,:),2^spin_num,2^spin_num);
    end

    U = reshape(ULs(2*layer,:,:),2^spin_num,2^spin_num);

    URs = complex(zeros([2*layer,2^spin_num,2^spin_num],'single'));
    URs(2*layer,:,:) = cz*reshape(Ry(layer,:,:),2^spin_num,2^spin_num);
    URs(2*layer-1,:,:) = reshape(URs(2*layer,:,:),2^spin_num,2^spin_num)*reshape(Rx(layer,:,:),2^spin_num,2^spin_num);
%     U = eye(2^spin_num,'single');
    for i = linspace(layer-1,1,layer-1)
        URs(2*i,:,:) = reshape(URs(2*i+1,:,:),2^spin_num,2^spin_num)*cz*reshape(Ry(i,:,:),2^spin_num,2^spin_num);
        URs(2*i-1,:,:) = reshape(URs(2*i,:,:),2^spin_num,2^spin_num)*reshape(Rx(i,:,:),2^spin_num,2^spin_num);
    end
    
    U_diff_b = complex(zeros([2*spin_num*layer,2^spin_num,2^spin_num],'single'));
    for i = 1:layer
        if i == 1
          for j =1:spin_num
            U_diff_b(2*spin_num*(i-1)+j,:,:) = -0.5*1j*reshape(URs(i,:,:),2^spin_num,2^spin_num)*sigmax_multi(spin_num,j);
            U_diff_b(2*spin_num*(i-1)+spin_num+j,:,:) = -0.5*1j*reshape(URs(2*i,:,:),2^spin_num,2^spin_num)*sigmay_multi(spin_num,j)*reshape(ULs(2*i-1,:,:),2^spin_num,2^spin_num);
          end
        else
            for j = 1:spin_num
                U_diff_b(2*spin_num*(i-1)+j,:,:) = -0.5*1j*reshape(URs(2*i-1,:,:),2^spin_num,2^spin_num)*sigmax_multi(spin_num,j)*reshape(ULs(2*i-2,:,:),2^spin_num,2^spin_num);
                U_diff_b(2*spin_num*(i-1)+spin_num+j,:,:) = -0.5*1j*reshape(URs(2*i,:,:),2^spin_num,2^spin_num)*sigmay_multi(spin_num,j)*reshape(ULs(2*i-1,:,:),2^spin_num,2^spin_num);
            end
        end
    end
end

function [U, U_diff_b] = generate_Trotter_ansatz(spin_num, layer, bvalue)
    Rxs = complex(zeros([spin_num*layer,2,2],'single'));
    Rzzs = complex(zeros([(spin_num-1)*layer,2^spin_num,2^spin_num],'single'));
    for i = 1:layer
        for j = 1:spin_num
            Rxs(spin_num*(i-1)+j,:,:) = rx(bvalue((2*spin_num-1)*(i-1)+j));
            if j<spin_num
                Rzzs((spin_num-1)*(i-1)+j,:,:) = rzz_multi(spin_num,j,bvalue((2*spin_num-1)*(i-1)+j+spin_num));
            end
        end
    end
    Rx = complex(zeros([layer,2^spin_num,2^spin_num],'single'));
    Rzz1 = complex(zeros([layer,2^spin_num,2^spin_num],'single'));
    Rzz2 = complex(zeros([layer,2^spin_num,2^spin_num],'single'));
    if spin_num>2
        for i = 1:layer
            x = reshape(Rxs(spin_num*(i-1)+1,:,:),2,2);
            zz1 = reshape(Rzzs((spin_num-1)*(i-1)+1,:,:),2^spin_num,2^spin_num);
            zz2 = reshape(Rzzs((spin_num-1)*(i-1)+2,:,:),2^spin_num,2^spin_num);
            for j = 2:spin_num
                x = kron(x,reshape(Rxs(spin_num*(i-1)+j,:,:),2,2));
            end
            Rx(i,:,:) = x;
            for j = 3:spin_num-1
                if mod(j,2) == 1
                    zz1 = reshape(Rzzs((spin_num-1)*(i-1)+j,:,:),2^spin_num,2^spin_num)*zz1; 
                else
                    zz2 = reshape(Rzzs((spin_num-1)*(i-1)+j,:,:),2^spin_num,2^spin_num)*zz2; 
                end
            end
            Rzz1(i,:,:) = zz1;
            Rzz2(i,:,:) = zz2;
        end
        ULs = complex(zeros([3*layer,2^spin_num,2^spin_num],'single'));
        ULs(1,:,:) = reshape(Rx(1,:,:),2^spin_num,2^spin_num);
        ULs(2,:,:) = reshape(Rzz1(1,:,:),2^spin_num,2^spin_num)*reshape(ULs(1,:,:),2^spin_num,2^spin_num);
        ULs(3,:,:) = reshape(Rzz2(1,:,:),2^spin_num,2^spin_num)*reshape(ULs(2,:,:),2^spin_num,2^spin_num);
        for i = 2:layer
            ULs(3*i-2,:,:) = reshape(Rx(i,:,:),2^spin_num,2^spin_num)*reshape(ULs(3*i-3,:,:),2^spin_num,2^spin_num);
            ULs(3*i-1,:,:) = reshape(Rzz1(i,:,:),2^spin_num,2^spin_num)*reshape(ULs(3*i-2,:,:),2^spin_num,2^spin_num);
            ULs(3*i,:,:) = reshape(Rzz2(i,:,:),2^spin_num,2^spin_num)*reshape(ULs(3*i-1,:,:),2^spin_num,2^spin_num);
        end

        U = reshape(ULs(3*layer,:,:),2^spin_num,2^spin_num);

        URs = complex(zeros([3*layer,2^spin_num,2^spin_num],'single'));
        URs(3*layer,:,:) = reshape(Rzz2(layer,:,:),2^spin_num,2^spin_num);
        URs(3*layer-1,:,:) = reshape(URs(3*layer,:,:),2^spin_num,2^spin_num)*reshape(Rzz1(layer,:,:),2^spin_num,2^spin_num);        
        URs(3*layer-2,:,:) = reshape(URs(3*layer-1,:,:),2^spin_num,2^spin_num)*reshape(Rx(layer,:,:),2^spin_num,2^spin_num);
        for i = linspace(layer-1,1,layer-1)
            URs(3*i,:,:) = reshape(URs(3*i+1,:,:),2^spin_num,2^spin_num)*reshape(Rzz2(i,:,:),2^spin_num,2^spin_num);
            URs(3*i-1,:,:) = reshape(URs(3*i,:,:),2^spin_num,2^spin_num)*reshape(Rzz1(i,:,:),2^spin_num,2^spin_num);
            URs(3*i-2,:,:) = reshape(URs(3*i-1,:,:),2^spin_num,2^spin_num)*reshape(Rx(i,:,:),2^spin_num,2^spin_num);
        end
    
        U_diff_b = complex(zeros([(2*spin_num-1)*layer,2^spin_num,2^spin_num],'single'));
        for i = 1:layer
            if i == 1
                for j =1:spin_num
                    U_diff_b((2*spin_num-1)*(i-1)+j,:,:) = -0.5*1j*reshape(URs(i,:,:),2^spin_num,2^spin_num)*sigmax_multi(spin_num,j);
                    if j<spin_num
                        if mod(j,1) == 1
                            U_diff_b((2*spin_num-1)*(i-1)+spin_num+j,:,:) = -0.5*1j*reshape(URs(3*i-1,:,:),2^spin_num,2^spin_num)*sigmaz_sigmaz_multi(spin_num,j)*reshape(ULs(3*i-2,:,:),2^spin_num,2^spin_num);
                        else
                             U_diff_b((2*spin_num-1)*(i-1)+spin_num+j,:,:) = -0.5*1j*reshape(URs(3*i,:,:),2^spin_num,2^spin_num)*sigmaz_sigmaz_multi(spin_num,j)*reshape(ULs(3*i-1,:,:),2^spin_num,2^spin_num);
                        end
                    end
                end
            else
                for j = 1:spin_num
                    U_diff_b((2*spin_num-1)*(i-1)+j,:,:) = -0.5*1j*reshape(URs(3*i-2,:,:),2^spin_num,2^spin_num)*sigmax_multi(spin_num,j)*reshape(ULs(3*i-3,:,:),2^spin_num,2^spin_num);
                    if j<spin_num
                        if mod(j,1) == 1
                            U_diff_b((2*spin_num-1)*(i-1)+spin_num+j,:,:) = -0.5*1j*reshape(URs(3*i-1,:,:),2^spin_num,2^spin_num)*sigmaz_sigmaz_multi(spin_num,j)*reshape(ULs(3*i-2,:,:),2^spin_num,2^spin_num);
                        else
                             U_diff_b((2*spin_num-1)*(i-1)+spin_num+j,:,:) = -0.5*1j*reshape(URs(3*i,:,:),2^spin_num,2^spin_num)*sigmaz_sigmaz_multi(spin_num,j)*reshape(ULs(3*i-1,:,:),2^spin_num,2^spin_num);
                        end
                    end
                 end
            end
        end
    end
    if spin_num == 2
        for i = 1:layer
            x = reshape(Rxs(spin_num*(i-1)+1,:,:),2,2);
            zz1 = reshape(Rzzs(i,:,:),2^spin_num,2^spin_num);
            for j = 2:spin_num
                x = kron(x,reshape(Rxs(spin_num*(i-1)+j,:,:),2,2));
            end
            Rx(i,:,:) = x;
            Rzz1(i,:,:) = zz1;
        end
        ULs = complex(zeros([2*layer,2^spin_num,2^spin_num],'single'));
        ULs(1,:,:) = reshape(Rx(1,:,:),2^spin_num,2^spin_num);
        ULs(2,:,:) = reshape(Rzz1(1,:,:),2^spin_num,2^spin_num)*reshape(ULs(1,:,:),2^spin_num,2^spin_num);
        for i = 2:layer
            ULs(2*i-1,:,:) = reshape(Rx(i,:,:),2^spin_num,2^spin_num)*reshape(ULs(2*i-2,:,:),2^spin_num,2^spin_num);
            ULs(2*i,:,:) = reshape(Rzz1(i,:,:),2^spin_num,2^spin_num)*reshape(ULs(2*i-1,:,:),2^spin_num,2^spin_num);
        end

        U = reshape(ULs(2*layer,:,:),2^spin_num,2^spin_num);

        URs = complex(zeros([2*layer,2^spin_num,2^spin_num],'single'));
        URs(2*layer,:,:) = reshape(Rzz1(layer,:,:),2^spin_num,2^spin_num);
        URs(2*layer-1,:,:) = reshape(URs(2*layer,:,:),2^spin_num,2^spin_num)*reshape(Rx(layer,:,:),2^spin_num,2^spin_num);        
        for i = linspace(layer-1,1,layer-1)
            URs(2*i,:,:) = reshape(URs(2*i+1,:,:),2^spin_num,2^spin_num)*reshape(Rzz1(i,:,:),2^spin_num,2^spin_num);
            URs(2*i-1,:,:) = reshape(URs(2*i,:,:),2^spin_num,2^spin_num)*reshape(Rx(i,:,:),2^spin_num,2^spin_num);
        end
    
        U_diff_b = complex(zeros([(2*spin_num-1)*layer,2^spin_num,2^spin_num],'single'));
        for i = 1:layer
            if i == 1
                for j =1:spin_num
                    U_diff_b((2*spin_num-1)*(i-1)+j,:,:) = -0.5*1j*reshape(URs(1,:,:),2^spin_num,2^spin_num)*sigmax_multi(spin_num,j);  
                end
                U_diff_b((2*spin_num-1)*(i-1)+spin_num+1,:,:) = -0.5*1j*reshape(URs(2*i,:,:),2^spin_num,2^spin_num)*sigmaz_sigmaz_multi(spin_num,1)*reshape(ULs(2*i-1,:,:),2^spin_num,2^spin_num);
            else
                for j = 1:spin_num
                    U_diff_b((2*spin_num-1)*(i-1)+j,:,:) = -0.5*1j*reshape(URs(2*i-1,:,:),2^spin_num,2^spin_num)*sigmax_multi(spin_num,j)*reshape(ULs(2*i-2,:,:),2^spin_num,2^spin_num);
                end
                 U_diff_b((2*spin_num-1)*(i-1)+spin_num+1,:,:) = -0.5*1j*reshape(URs(2*i,:,:),2^spin_num,2^spin_num)*sigmaz_sigmaz_multi(spin_num,1)*reshape(ULs(2*i-1,:,:),2^spin_num,2^spin_num);
            end
        end
    end
end

function phi = computational_basis(spin_num,index)
    phi = zeros([2^spin_num,1],'single');
    phi(index) = 1;
end

function [rho, rho_diff_a, rho_diff_b] = rho_state(U, U_diff_b, R, spin_num, avalue, bvalue)
    num_a = R;
    num_b = size(bvalue);
    num_b = num_b(2);
    Utotal = kron(U,conj(U));
    phi_a = complex(zeros([2^(2*spin_num),1],'single'));
    for i = 1:R
        phi_a = phi_a + avalue(i)*kron(computational_basis(spin_num,i),computational_basis(spin_num,i));
    end
    rho = Utotal*phi_a;
    rho_diff_a = complex(zeros([2^(2*spin_num),num_a],'single'));
    rho_diff_b = complex(zeros([2^(2*spin_num),num_b],'single'));
    for i = 1:num_a
        rho_diff_a(:,i) = Utotal*kron(computational_basis(spin_num,i),computational_basis(spin_num,i));
    end
    for i = 1:num_b
        U_diff = reshape(U_diff_b(i,:,:),2^spin_num,2^spin_num);
        rho_diff_b(:,i) = (kron(U_diff,conj(U))+kron(U,conj(U_diff)))*phi_a;
    end
end

function L = generate_Liou_TFIM(spin_num, Jz, h, gamma)
    sigmaz = [[1,0];[0,-1]];
    sigmax = [[0,1];[1,0]];
    sigmam = [[0,0];[1,0]];

    Hzz = zeros([spin_num,2^spin_num,2^spin_num],'single');
    zz = kron(sigmaz,sigmaz);
    for i = 1:spin_num
        if i == spin_num
            hzz = kron(sigmaz,eye(2^(spin_num-2),'single'));
            hzz = kron(hzz,sigmaz);
            Hzz(i,:,:) = hzz;
        else
            hzz = kron(eye(2^(i-1),'single'),zz);
            hzz = kron(hzz,eye(2^(spin_num-i-1),'single'));
            Hzz(i,:,:) = hzz;
        end
    end

    Hx = zeros([spin_num,2^spin_num,2^spin_num],'single');
    for i = 1:spin_num
        hx = kron(eye(2^(i-1),'single'),sigmax);
        hx = kron(hx,eye(2^(spin_num-i),'single'));        
        Hx(i,:,:) = hx;
    end

    H = 0;
    if spin_num == 2
        H = H+Jz*Hzz(1,:,:)+h*Hx(1,:,:)+h*Hx(2,:,:);
    else
        for i = 1:spin_num
            H = H+Jz*Hzz(i,:,:)+h*Hx(i,:,:);
        end
        %H = H-Jz*Hzz(spin_num,:,:);
    end
    H = reshape(H,2^spin_num,2^spin_num);
    cops = zeros([spin_num,2^spin_num,2^spin_num],'single');
    for i = 1:spin_num
        cop = kron(eye(2^(i-1),'single'),sigmam);
        cop = kron(cop,eye(2^(spin_num-i),'single'));               
        cops(i,:,:) = cop;
    end    
    L = -1j*(kron(H,eye(2^spin_num,'single'))-kron(eye(2^spin_num,'single'),transpose(H)));
    for i = 1:spin_num
        cop = reshape(cops(i,:,:),2^spin_num,2^spin_num);
        L = L+gamma*(kron(cop,cop)-0.5*kron(cop'*cop,eye(2^spin_num,'single'))-0.5*kron(eye(2^spin_num,'single'),cop'*cop));
    end
end

function M = matrix_elements(rho_diff_a, rho_diff_b, R)
    num_a = R;
    num_b = size(rho_diff_b);
    num_b = num_b(2);
    num = num_a + num_b;
    M = zeros([num,num],'single');
    for k = 1:num_a
        for j = 1:num_a
            if  k == j 
                M(k,j) = 1;
            else 
                M(k,j) = 0;
            end
        end

        for j = 1:num_b
            M(k,num_a+j) = real(rho_diff_a(:,k)'*rho_diff_b(:,j));
            M(num_a+j,k) =  M(k,num_a+j);
        end
    end

    for k = 1:num_b
        for j = 1:num_b
            M(num_a+k,num_a+j) = real(rho_diff_b(:,k)'*rho_diff_b(:,j));
        end
    end
    M = single(M);
end

function V = vec_elements(rho, rho_diff_a, rho_diff_b, L, R)
    num_a = R;
    num_b = size(rho_diff_b);
    num_b = num_b(2);
    num = num_a + num_b;
    V = zeros([num,1],'single');
    for i = 1:num_a
        V(i) = real(rho_diff_a(:,i)'*L*rho);
    end
    for i = 1:num_b
        V(num_a+i) = real(rho_diff_b(:,i)'*L*rho);
    end
    V = single(V);
end

function mz = Average_Magnetization(U,spin_num,R,avalue)
    rho_eig = avalue/sum(avalue);
    num_a = R;
    Oz = 0; 
    for i = 1:spin_num
        Oz  = Oz+sigmaz_multi(spin_num,i)/spin_num;
    end
    Oz = reshape(Oz,2^spin_num,2^spin_num);
    mz = 0;
    for i = 1:num_a
        phi = computational_basis(spin_num,i);
        phi = U*phi;
        mz = mz+ (phi'*Oz*phi)*rho_eig(i);
    end
    mz = single(mz);
end

function mx = Average_Magnetization_X(U,spin_num,R,avalue)
    rho_eig = avalue/sum(avalue);
    num_a = R;
    Ox = 0; 
    for i = 1:spin_num
        Ox  = Ox+sigmax_multi(spin_num,i)/spin_num;
    end
    Ox = reshape(Ox,2^spin_num,2^spin_num);
    mx = 0;
    for i = 1:num_a
        phi = computational_basis(spin_num,i);
        phi = U*phi;
        mx = mx+ (phi'*Ox*phi)*rho_eig(i);
    end
    mx = single(mx);
end

function diff = diff_para(rho_value, rho_diff_a_value, rho_diff_b_value, L, R)
    M = matrix_elements(rho_diff_a_value, rho_diff_b_value, R);
    V = vec_elements(rho_value, rho_diff_a_value, rho_diff_b_value, L, R);
    diff = lsqr(M,V,10^(-6),100);
    %diff = linsolve(M,V);
end

function [tlist, mz, mx, paras] = VQD(L, R, spin_num, para_initial, dt, steps, layer)
    mz = zeros([1,steps],'single');
    mx = zeros([1,steps],'single');
    tlist = zeros([1,steps],'single');
    num_b = size(para_initial);
    num_b = num_b(2)-R;
    paras = zeros([R+num_b,steps],'single');
    paras(:,1) = para_initial';
    for i = 1:steps
        %[U, U_diff_b] = generate_hardware_efficient_ansatz(spin_num, layer, paras(R+1:R+num_b,i)');
        [U, U_diff_b] = generate_Trotter_ansatz(spin_num, layer, paras(R+1:R+num_b,i)');
        [rho, rho_diff_a, rho_diff_b] = rho_state(U, U_diff_b, R, spin_num, paras(1:R,i)', paras(R+1:R+num_b,i)');       
        tlist(i) = (i-1)*dt;
        mz(i) = Average_Magnetization(U,spin_num,R,paras(1:R,i)');
        mx(i) = Average_Magnetization_X(U,spin_num,R,paras(1:R,i)');
        disp(tlist(i));
        disp(mz(i));
        if i<steps
            diff = diff_para(rho, rho_diff_a, rho_diff_b, L, R);
            %disp(diff);
            paras(:,i+1) = paras(:,i)+diff*dt;
        end
    end
%     figure
%     plot(tlist,mz,'b-.');
%     saveas(gcf,'mz.jpg');
%     figure
%     plot(tlist,mx,'b-.');
%     saveas(gcf,'mx.jpg');
end 