function plot_state_vars(tarray,zarray)
t= tarray;
q1 = zarray(:,1); q2 = zarray(:,2); q3 = zarray(:,3); q4 = zarray(:,4);
w1 = zarray(:,5); w2 = zarray(:,6); w3 = zarray(:,7);

qarray = [q1';q2';q3';q4'];
b3array=zeros(3,1,length(qarray));
for i=1:length(qarray)
b3array(:,:,i) = rotateframe(qarray(:,i),[0;0;1]);
end
b3array

% Graph 1: Angular Velocities
figure(1)
clf
subplot(2,1,1)
plot(t,w1,'b--',t,w2,'r--')
set(gca,'FontName','Times','FontSize',12)
legend({'\omega_1','\omega_2'})
ylabel('\omega_1, \omega_2 (rad/s)')
xlim([0,t(end)])
subplot(2,1,2)
plot(t,w3,'--')
set(gca,'FontName','Times','FontSize',12)
ylabel('\omega_3 (rad/s)')
ylim([w3(1)*0.99,w3(1)*1.01])
xlabel('Time (s)')
xlim([0,t(end)])

% Graph 2: Quaternions
figure()
subplot(2,1,1)
plot(t,q1,t,q2,t,q3,t,q4)
legend('q1','q2','q3','q4')

subplot(2,1,2)
quiver3(zeros(1,1,length(qarray)),zeros(1,1,length(qarray)),...
    zeros(1,1,length(qarray)),...
    b3array(1,:,:),b3array(2,:,:),b3array(3,:,:))
end