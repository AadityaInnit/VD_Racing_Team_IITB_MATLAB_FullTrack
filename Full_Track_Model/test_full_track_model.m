function tests = test_full_track_model()
%TEST_FULL_TRACK_MODEL Physical balances, planning, tracks and failure reporting.
tests = functiontests(localfunctions);
end

function setupOnce(t)
c = full_track_config(); c.makePlots=false;
data = load(fullfile(fileparts(mfilename('fullpath')),'data','four_wheel_ggv.mat'),'ggv');
t.TestData.c=c; t.TestData.L=data.ggv;
a=linspace(0,2*pi,81).';
t.TestData.ellipse=make_full_track([35*sin(a),22*(1-cos(a))],true,400);
t.TestData.out=full_track_model(c,t.TestData.ellipse,GGV=data.ggv);
end

function testRWDOpenDifferentialAndPower(t)
c=t.TestData.c;
d=probe(c,[20;0.5;0.4;0.08],1e6);
verifyEqual(t,d.Fx(1:2),[0;0],'AbsTol',1e-12);
verifyEqual(t,d.Fx(3),d.Fx(4),'AbsTol',1e-10);
verifyLessThanOrEqual(t,sum(d.Fx.*d.wheelLongitudinalSpeed),c.powerLimit+1e-6);
verifyLessThanOrEqual(t,d.utilization,1+1e-10);
end

function testIndividualWheelKinematicsAndYaw(t)
c=t.TestData.c; p=c.vehicle;
d=probe(c,[12;0.4;0.5;0.10],800);
verifyGreaterThan(t,d.wheelSteer(1),d.wheelSteer(2));
verifyGreaterThan(t,abs(d.alpha(1)-d.alpha(2)),1e-4);
verifyGreaterThan(t,abs(d.alpha(3)-d.alpha(4)),1e-4);
x=[p.a;p.a;-p.b;-p.b]; y=[p.track_f;-p.track_f;p.track_r;-p.track_r]/2;
verifyEqual(t,d.Mz,sum(x.*d.yBody-y.*d.xBody),'AbsTol',1e-10);
end

function testVerticalLoadAndLateralTransfer(t)
c=t.TestData.c; p=c.vehicle;
d=probe(c,[12;0;0.5;0.12],500);
verifyTrue(t,d.loadConverged);
verifyGreaterThan(t,d.FY,0);
verifyGreaterThan(t,d.Fz(2),d.Fz(1));
verifyGreaterThan(t,d.Fz(4),d.Fz(3));
verifyEqual(t,sum(d.Fz),p.m*p.g+0.5*p.rho*p.ClA*12^2,'AbsTol',1e-6);
verifyGreaterThanOrEqual(t,min(d.Fz),0);
end

function testYawInertiaChangesOnlyYawAcceleration(t)
c=t.TestData.c;
a=probe(c,[10;0.3;0.4;0.05],700); c.vehicle.Iz=2*c.vehicle.Iz;
b=probe(c,[10;0.3;0.4;0.05],700);
verifyEqual(t,a.Mz,b.Mz,'AbsTol',1e-10);
verifyEqual(t,a.Mz/t.TestData.c.vehicle.Iz,2*b.Mz/c.vehicle.Iz,'AbsTol',1e-10);
end

function testMatchedGGVBalances(t)
c=t.TestData.c; L=t.TestData.L; p=c.vehicle;
for mode=[1 -1]
    if mode>0, A=L.ax_acc; Z=L.state_acc; else, A=L.ax_brk; Z=L.state_brk; end
    [ia,iv]=find(isfinite(A));
    chosen=unique(round(linspace(1,numel(ia),min(15,numel(ia)))));
    for j=chosen
        z=reshape(Z(ia(j),iv(j),:),4,1); V=L.v(iv(j)); beta=z(1); ay=L.ay(ia(j),iv(j));
        d=probe(c,[V*cos(beta);V*sin(beta);ay/V;z(2)],z(3));
        at=(d.FX*cos(beta)+d.FY*sin(beta))/p.m;
        an=(-d.FX*sin(beta)+d.FY*cos(beta))/p.m;
        verifyLessThan(t,abs(at-z(4))/p.g,1e-4);
        verifyLessThan(t,abs(an-ay)/p.g,1e-4);
        verifyLessThan(t,abs(d.Mz)/(p.m*p.g*p.wb),1e-4);
    end
end
end

function testEllipseCompletesWithDynamicStates(t)
o=t.TestData.out;
verifyTrue(t,o.completed); verifyTrue(t,isfinite(o.lapTime));
verifyEqual(t,o.speed(1),0,'AbsTol',1e-12);
verifyEqual(t,o.s(end),o.track.length,'AbsTol',1e-6);
verifyGreaterThan(t,max(abs(o.Vy)),0.02);
verifyGreaterThan(t,max(abs(o.rDot)),0.02);
verifyLessThan(t,max(abs(o.e_y)),o.config.trackHalfWidth);
verifyLessThan(t,max(o.loadResidual),1e-3);
verifyLessThanOrEqual(t,max(o.tireUtilization),1+1e-10);
verifyGreaterThan(t,max(o.referenceCurvature)-min(o.referenceCurvature),0.02);
end

function testNewtonEquationsAndWorldTrajectory(t)
o=t.TestData.out; p=o.config.vehicle;
% Integrate Newton's equations independently; finite differences amplify
% force-allocation kinks at braking/traction saturation.
u=o.Vx(1)+cumtrapz(o.t,o.FX/p.m+o.r.*o.Vy);
v=o.Vy(1)+cumtrapz(o.t,o.FY/p.m-o.r.*o.Vx);
r=o.r(1)+cumtrapz(o.t,o.Mz/p.Iz);
verifyLessThan(t,max(abs(u-o.Vx)),0.02);
verifyLessThan(t,max(abs(v-o.Vy)),0.02);
verifyLessThan(t,max(abs(r-o.r)),0.01);
vx=o.Vx.*cos(o.psi)-o.Vy.*sin(o.psi); vy=o.Vx.*sin(o.psi)+o.Vy.*cos(o.psi);
X=o.X(1)+cumtrapz(o.t,vx); Y=o.Y(1)+cumtrapz(o.t,vy);
verifyLessThan(t,max(abs(X-o.X)),0.03);
verifyLessThan(t,max(abs(Y-o.Y)),0.03);
end

function testClockwiseFlyingLapAndSeam(t)
c=t.TestData.c; c.startMode="flying";
a=linspace(0,2*pi,81).'; track=make_full_track([-25*sin(a),25*(1-cos(a))],true,300);
o=full_track_model(c,track,GGV=t.TestData.L);
verifyTrue(t,o.completed);
verifyLessThan(t,mean(o.r),0);
verifyLessThan(t,o.plan.cyclicResidual,1e-5);
verifyLessThan(t,max(abs(o.e_y)),0.08);
end

function testOpenCourseWithBothTurnDirections(t)
c=t.TestData.c; c.speedScale=0.6;
x=linspace(0,160,81).'; track=make_full_track([x,5*sin(2*pi*x/100)],false,400);
verifyGreaterThan(t,max(track.curvature),0.01);
verifyLessThan(t,min(track.curvature),-0.01);
o=full_track_model(c,track,GGV=t.TestData.L);
verifyTrue(t,o.completed);
verifyGreaterThan(t,max(o.r),0.02); verifyLessThan(t,min(o.r),-0.02);
verifyLessThan(t,max(abs(o.e_y)),0.5);
end

function testRejectStaleGGVAndInvalidTrack(t)
c=t.TestData.c; c.frontRollTransferFraction=0.6;
verifyError(t,@()full_track_model(c,t.TestData.ellipse,GGV=t.TestData.L),'full_track:GGVMismatch');
verifyError(t,@()make_full_track(zeros(8,2),true),'full_track:Track');
end

function testFailureDoesNotProduceLapTime(t)
c=t.TestData.c; c.trackHalfWidth=0.001;
o=full_track_model(c,t.TestData.ellipse,GGV=t.TestData.L);
verifyFalse(t,o.completed); verifyTrue(t,isnan(o.lapTime));
verifyEqual(t,o.status,"failed: track boundary exceeded");
end

function testTimeStepConvergence(t)
c=t.TestData.c; c.maxStep=0.01;
o=full_track_model(c,t.TestData.ellipse,GGV=t.TestData.L);
verifyTrue(t,o.completed);
verifyLessThan(t,abs(o.lapTime-t.TestData.out.lapTime),0.005);
end

function d=probe(c,z,F)
d=full_track_model(c,[],PlantProbe=struct('state',z,'force',F));
end
