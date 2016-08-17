#include <pngwriter.h>

void plot_mesh(vector<cell*>& domain)
{
	pngwriter domain(300, 300, 0,"domain.png");

	for (int i = 0; i < domain.size(); i++)
		for (int j= 0; j< domain[i]->edges.size(); ++j)
			png.line(domain[i]->edges[j].s.first, domain[i]->edges[j].s.second, domain[i]->edges[j].e.first, domain[i]->edges[j].e.second, 1.0, 1.0, 1.0);
}
