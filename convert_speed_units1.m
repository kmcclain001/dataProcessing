
%% convert speed to cm/s for session 
% (little bit of an art here)
    
function behavior = convert_speed_units1(varargin)

% Parse Inputs

p = inputParser;
addParameter(p,'basepath',pwd,@isstr);
addParameter(p,'filename',[],@isstr);
addParameter(p,'trackLength',[],@isfloat);

parse(p,varargin{:});

basepath = p.Results.basepath;
filename = p.Results.filename;
track_length = p.Results.trackLength;

load([basepath filesep filename '.linear.behavior.mat']);

fps = behavior.samplingRate;

if isfield(behavior,'units')
    u = behavior.units;

    switch u
        case 'pixels'
            %cm_p_b = 1/1.6666;
            cm_p_b = 1/3;
            %count1 = count1 + 1;
        case 'm'
            cm_p_b = 100;
            %count2 = count2 + 1;
        case 'mm'
            cm_p_b = 1/10;
            %count3 = count3 + 1;
        otherwise
            error('some other spatial unit')
    end
    
    for i =1:length(behavior.events.trials)
        s = behavior.events.trials{i}.speedArb';
        s_new = s*cm_p_b*fps;
        behavior.events.trials{i}.speed = s_new';
    end
elseif ~isempty(track_length)
    x = behavior.position.x;
    y = behavior.position.y;
    
    maxx = max(x);
    minx = min(x);
    maxy = max(y);
    miny = min(y);
    
    fps = behavior.samplingRate;
    dx = maxx - minx;
    dy = maxy-miny;
    d = (dx.^2+dy.^2).^.5;
    
    for i = 1:length(behavior.events.trials)
        s = behavior.events.trials{i}.speedArb';
        s_new = s *(track_length/d)*fps;
        behavior.events.trials{i}.speed = s_new';
    end
    behavior.units = track_length/d;
else
    error('no units or tracklength known')
end

    save([basepath filesep filename '.linear.behavior.mat'],'behavior')
    
end
