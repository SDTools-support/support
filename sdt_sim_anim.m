function   out=sdt_sim_anim(t,u)

persistent cf 

if ischar(t)
[CAM,Cam]=comstr(t,1); 
if comstr(Cam,'init')
%% #Init commands
    
  ref=evalin('base','ref');
  if isfield(ref,'AnSampleTime') % dt,T0 set in Blockparameters
   r1=ref.AnSampleTime; if ~ischar(r1);r1=comstr(r1,-30);end
   set_param('sim_anim/anim_struct','SystemSampleTime',r1)
  end
  if isfield(ref,'cf')
   cf=ref.cf; % Actually declare the proper feplot
   figure(ref.cf.opt(1)); % Raise animation so that it can be viewed
  else; cf=feplot;
  end
  cf.data.LastWhen=now; 
   
else
%%
    error('%s unkown',CAM)
end
    
elseif nargin==2 
%% #anim t>0 % Actually do an animation
   if isempty(cf);cf=feplot; end    
   if ~isa(cf,'sdth')
   elseif strcmp(get(cf.opt(1),'tag'),'feplot')
    %cf=get(2,'userdata');
    if now-cf.data.LastWhen*24*3600>.1 % Only animate every .1 s
     of_time(-1,cf.def.def,u,zeros(1));
     cf.def.data(1)=t; feplot(cf);
    end
   end
   
else % Do something else
  out=zeros(4,1);
end
