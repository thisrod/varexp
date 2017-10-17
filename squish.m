% basic squeezed state play

% The small multiples form a half plane: squeezing in the opposite direction has the same effect.

global N; N = 15; brackets
mesh = -5:0.1:5;  phasespace

figure, zplot(mesh, mesh, Aps(end,:)), axis image
title 'Maximum number state'

figure, zplot(mesh, mesh, Aps'*nq*evan(0, 'even')), axis image
title 'Numerically exact vacuum'

figure, zplot(mesh, mesh, Aps'*sqn(0,0.1)), axis image
title 'Minimally squeezed vacuum'

figure, zplot(mesh, mesh, Aps'*sqn(0,1)), axis image
title 'Squeezed vacuum'

figure, zplot(mesh, mesh, Aps'*sqn(0,0.7i)), axis image
title 'Squeezed vacuum, +i direction'

figure, zplot(mesh, mesh, Aps'*sqn(0,-0.7i)), axis image
title 'Squeezed vacuum, -i direction'

figure, zplot(mesh, mesh, Aps'*sqn(0,-0.5-0.7i)), axis image
title 'Squeezed vacuum, -1-i direction'

figure, zplot(mesh, mesh, Aps'*sqn(0,10)), axis image
title 'Very squeezed vacuum'

figure, zplot(mesh, mesh, Aps'*sqn(exp(1)*2,1)), axis image
title 'Squeezed state with coherent amplitude 2'
