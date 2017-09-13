function show_coins(centroids,types,total_value)

% show_coins
%
% centroids: (k x 2) array of (x,y) centroids
% types: k x 1 array of coin type (1: 50 cents, 2: 20 cents, 3: 10 cents, 4: 5 cents)
% total_value: total value in dollars
%

labels = {'50','20','10','5'};

hold on
for j = 1:length(types)
    h1 = plot(centroids(j,1),centroids(j,2),'r.');
    set(h1,'MarkerSize',15)
    t1 = text(centroids(j,1),centroids(j,2),labels{types(j)});
    set(t1,'FontSize',16);
    set(t1,'Color',[0,1,0]);
end

t_total = text(10,10,['coins: ',num2str(length(types)),' value: $ ',num2str(total_value)]);
set(t_total,'FontSize',16);
set(t_total,'Color',[1,0,0]);

hold off
