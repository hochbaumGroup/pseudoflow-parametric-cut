from abc import ABC, abstractmethod


class GraphWrapper(ABC):
    def __init__(self, G):
        self.G = G

    @abstractmethod
    def predecessors(self, i):
        pass

    @abstractmethod
    def successors(self, i):
        pass

    @abstractmethod
    def edge_attribute_value(self, i, j, attribute):
        pass

    @abstractmethod
    def num_nodes(self):
        pass

    @abstractmethod
    def num_edges(self):
        pass

    @abstractmethod
    def edges(self, data=False):
        pass

    @abstractmethod
    def nodes(self):
        pass


class NetworkxGraphWrapper(GraphWrapper):
    def nodes(self):
        return self.G.nodes()

    def predecessors(self, i):
        return self.G.predecessors(i)

    def edge_attribute_value(self, i, j, attribute):
        return self.G[i][j][attribute]

    def num_nodes(self):
        return self.G.number_of_nodes()

    def num_edges(self):
        return self.G.number_of_edges()

    def edges(self, data=False):
        return self.G.edges(data=data)

    def successors(self, i):
        return self.G.successors(i)


class IgraphGraphWrapper(GraphWrapper):
    def nodes(self):
        for v in self.G.vs:
            yield v.index

    def predecessors(self, i):
        return self.G.predecessors(i)

    def successors(self, i):
        return self.G.successors(i)

    def edge_attribute_value(self, i, j, attribute):
        return self.G.es[self.G.get_eid(i, j)][attribute]

    def num_nodes(self):
        return self.G.vcount()

    def num_edges(self):
        return self.G.ecount()

    def edges(self, data=False):
        if data:
            for e in self.G.es:
                yield e.source, e.target, e
        else:
            for e in self.G.es:
                yield e.source, e.target
