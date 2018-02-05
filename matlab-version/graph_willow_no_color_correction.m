function graph_willow_no_color_correction(m,z,q,t,h,T,G,P)

% GRAPHICAL REPRESENTATION OF THE WILLOW TREE

% [A] BUILD GRID W

initial = [zeros(m,1) (t(2))*ones(m,1) G(:,1) G(:,2) q ones(m,1) zeros(m,2)];

start_date=kron(t(2:end-1), ones(m^2,1)); start_date=start_date(:);
end_date=kron(t(3:end), ones(m^2,1));  end_date=end_date(:);
start_node=kron(G(:,2:end-1), ones(m,1)); start_node=start_node(:);
end_node=kron(G(:,3:end), ones(1,m)); end_node=end_node(:);
transition=cell2mat(P);
transition=transition(:);

red=kron(rand(floor(m/2),numel(h)), ones(m,1));
green=kron(rand(floor(m/2),numel(h)), ones(m,1)); 
blue=kron(rand(floor(m/2),numel(h)), ones(m,1)); 

if mod(m, 2) == 0
red = [red; flip(red);]; red=red(:);
green = [green; flip(green)]; green=green(:);
blue = [blue; flip(blue)]; blue=blue(:);
else
red_mid = kron(rand(1,numel(h)), ones(m,1));
green_mid = kron(rand(1,numel(h)), ones(m,1)); 
blue_mid = kron(rand(1,numel(h)), ones(m,1)); 
red = [red; red_mid; flip(red)]; red=red(:);
green = [green; green_mid; flip(green)]; green=green(:);
blue = [blue; blue_mid; flip(blue)]; blue=blue(:);
end

W = [initial;start_date end_date start_node end_node transition red green blue];

% [B] PLOT THE RESULTS

steps=t*360;
steps_100 = linspace(0,T,200);
square_root = sqrt(steps_100).*z(m);

figure  
rgbcolor = [0 0.447 0.741];
plot(t,G,'.k');
whitebg([0.9451 0.9686 0.9490]);

hold on
set(gca,'XTick',t,'XtickLabel',steps,'Ydir','reverse');
xlabel('Time (days)'); ylabel('Spatial dimension');
hold on
for i=1:length(W)
if W(i,5) > 0
plot([W(i,1) W(i,2)],[W(i,3) W(i,4)]);
end
end
hold on
plot(steps_100,square_root,steps_100,-square_root,'Color',rgbcolor);
hold off
toc
end

