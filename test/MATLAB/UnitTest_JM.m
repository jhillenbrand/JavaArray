%% array definition
ar1 = rand(10e6, 1);
ar2 = rand(10e6, 1);

%% test scalar prod
t1 = cputime;
s = dot(ar1, ar2);
t_el = cputime - t1;

disp(['elapsed [s]: ' num2str(t_el)]);

%%