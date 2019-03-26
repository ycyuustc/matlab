clear all;
clear classes;

local = Local(1,100,4);
% sub = Sub(1,2,'truncate');
sub = Sub(local);
disp(sub);
tm = TransMat(local,sub);
disp(tm);

% construct super block
super = Super(local,sub,sub,tm,tm,tm,tm);
disp(super);
tmsys = TransMat(local,super.SysNew,super.SysDown,super.SysUp,'sys');
disp(tmsys);
tmenv = TransMat(local,super.EnvNew,super.EnvDown,super.EnvUp,'env');
disp(tmenv);

% truncate the tmsys
sub_trun_up = Sub(local,tmsys.Sub_up,super,'sys','up');
disp(sub_trun_up);
sub_trun_down = Sub(local,tmsys.Sub_down,super,'sys','down');
disp(sub_trun_down);
tmsys = TransMat(local,tmsys,sub_trun_up,sub_trun_down);
disp(tmsys);

% truncate the tmenv
sub_trun_up = Sub(local,tmenv.Sub_up,super,'sys','up');
disp(sub_trun_up);
sub_trun_down = Sub(local,tmenv.Sub_down,super,'sys','down');
disp(sub_trun_down);
tmenv = TransMat(local,tmenv,sub_trun_up,sub_trun_down);
disp(tmenv);

super = Super(local,super.SysNew,super.EnvNew,tmsys,tmenv);
disp(super);

tmsysup = TransMat(local,super.SysNew,tmsys,'sysup');
disp(tmsysup);

tmsysdown = TransMat(local,super.SysNew,tmsys,'sysdown');
disp(tmsysdown);

tmenvup = TransMat(local,super.SysNew,tmenv,'envup');
disp(tmenvup);

tmenvdown = TransMat(local,super.SysNew,tmenv,'envdown');
disp(tmenvdown);

sub_trun_up = Sub(local,tmsysup.Sub_up,super,'sys','up');
tmsysup = TransMat(local,tmsysup,sub_trun_up,'up');
disp(tmsysup);

sub_trun_down = Sub(local,tmsysup.Sub_up,super,'sys','up');
tmsysup = TransMat(local,tmsysup,sub_trun_up,'up');
disp(tmsysup);

% disp(sub);
% sub = Sub(sub,local,'add');
% disp(sub);
% sub = Sub(sub,local,'add');
% disp(sub);
% sub = Sub(sub,local,'add');
% disp(sub);

% tm = TransMat(local,sub);
% subBig = Sub(sub,local,'add');
% disp(tm);
% sup = Super(local,sub,sub,tm,tm,tm,tm);
% disp(sup);
% tms = TransMat(local,subBig,tm,tm,'sys');
% disp(tms);
% tme = TransMat(local,subBig,tm,tm,'env');
% disp(tme);
% 
% sup = Super(local,subBig,subBig,tms,tme);
% disp(sup);

% subBigBig = Sub(subBig,local,'add');
% disp(subBigBig);
% 
% sysdown = TransMat(local,subBigBig,tms,'sysdown');
% disp(sysdown);
% 
% sysup = TransMat(local,subBigBig,tms,'sysup');
% disp(sysup);
% 
% envdown = TransMat(local,subBigBig,tme,'envdown');
% disp(envdown);
% 
% envup = TransMat(local,subBigBig,tme,'envup');
% disp(envup);


