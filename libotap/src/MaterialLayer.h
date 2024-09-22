#pragma once

#include "OTypes.h"
#include "Material.h"
#include "mlinterp/mlinterp.hpp"
#include <numeric>
#include <memory>

namespace OTAP
{
    class Node
    {
    private:
        std::shared_ptr<Material> m_Material;

    public:
        Node(std::shared_ptr<Material> material, double tinit = 300.0) : m_Material(material), T(tinit) {}
        ~Node() = default;
        Node(const Node &) = default;

        double T;
        double k() { return m_Material->k[T]; }
        double rho() { return m_Material->rho[T]; }
        double rhochar() { return m_Material->rhochar[T]; }
        double rhopyrogas() { return m_Material->rhopyrogas[T]; }
        double Cp() { return m_Material->Cp[T]; }
        double Cpchar() { return m_Material->Cpchar[T]; }
        double Cppyrogas() { return m_Material->Cppyrogas[T]; }
        double emissivity() { return m_Material->emissivity[T]; }
        double Tabl() { return m_Material->Tabl; }
        double Tpyro() { return m_Material->Tpyro; }
        double Habl() { return m_Material->Habl; }
        double Hpyro() { return m_Material->Hpyro; }
        double rhovirgin() { return rho(Tpyro()); }
        double k(double pT) { return m_Material->k[pT]; }
        double rho(double pT) { return m_Material->rho[pT]; }
        double rhochar(double pT) { return m_Material->rhochar[pT]; }
        double rhopyrogas(double pT) { return m_Material->rhopyrogas[pT]; }
        double Cp(double pT) { return m_Material->Cp[pT]; }
        double Cpchar(double pT) { return m_Material->Cpchar[pT]; }
        double Cppyrogas(double pT) { return m_Material->Cppyrogas[pT]; }
        double emissivity(double pT) { return m_Material->emissivity[pT]; }
    };

    struct Layer
    {
        std::shared_ptr<Material> material;
        double thickness;
        size_t numNodes;

        Layer(std::shared_ptr<Material> mat, double thickness, size_t numNodes = 10)
            : material(mat), thickness(thickness), numNodes(std::max(size_t(2), numNodes))
        {
        }
        ~Layer() = default;
    };

    class LayerStack
    {
    private:
        std::vector<Layer> m_Layers;

    public:
        LayerStack() = default;
        LayerStack(std::initializer_list<Layer> layers) : m_Layers(layers) {}
        void Push(const Layer &layer)
        {
            m_Layers.push_back(layer);
        }
        template <typename... Args>
        void Emplace(Args &&...args)
        {
            m_Layers.emplace_back(std::forward<Args>(args)...);
        }
        void Delete(const size_t index)
        {
            m_Layers.erase(m_Layers.begin() + index);
        }
        auto GetCount() const { return m_Layers.size(); }
        auto &GetComponents() const { return m_Layers; }
        Layer &operator[](size_t index) noexcept { return m_Layers[index]; }
        Layer &at(size_t index) { return m_Layers.at(index); }
    };

    /// @brief Material layout from outside to inside
    class LayerMesh
    {
    private:
        std::shared_ptr<LayerStack> m_layers;

    public:
        bool NodesAreValid = false;
        std::vector<Node> Nodes;
        std::vector<size_t> NumNodes;

        LayerMesh(std::shared_ptr<LayerStack> ls) : m_layers(ls) {}
        LayerMesh(LayerStack ls) : m_layers(new LayerStack(ls)) {}
        ~LayerMesh() = default;
        auto GetCount() const { return m_layers->GetCount(); }
        auto &GetComponents() const { return m_layers->GetComponents(); }
        Layer &operator[](size_t index) noexcept { return m_layers->operator[](index); }

    private:
        void InitTemperature(double Tinit)
        {
            if (!NodesAreValid)
            {
                CreateNodes();
            }
            for (auto &&i : Nodes)
            {
                i.T = Tinit;
            }
        }

        void InitTemperature(TimeSeries Tinit)
        {
            auto &m_Layers = *m_layers;
            if (!NodesAreValid)
            {
                CreateNodes();
            }
            assert(m_Layers.GetCount() > 0);
            assert(Tinit.size() == (m_Layers.GetCount() + 1));

            for (size_t j = NumNodes[0]; j < NumNodes[1]; j++)
            {
                Nodes[j].T = Tinit[0] + (Tinit[1] - Tinit[0]) * (j - NumNodes[0]) / (NumNodes[1] - 1 - NumNodes[0]);
            }

            for (size_t i = 1; i < m_Layers.GetCount(); i++)
            {
                for (size_t j = NumNodes[i]; j < NumNodes[i + 1]; j++)
                {
                    Nodes[j].T = Tinit[i] + (Tinit[i + 1] - Tinit[i]) * (j - NumNodes[i] + 1) / (NumNodes[i + 1] - NumNodes[i]);
                }
            }
        }

    public:
        void CreateNodes()
        {
            // First Layer (one additional node)
            auto &m_Layers = *m_layers;
            assert(m_Layers.GetCount() > 0);
            for (size_t j = 0; j < m_Layers[0].numNodes + 1; j++)
            {
                Nodes.emplace_back(m_Layers[0].material);
            }
            NumNodes.push_back(0);
            NumNodes.push_back(m_Layers[0].numNodes + 1);

            // Subsequent Layer
            for (size_t i = 1; i < m_Layers.GetCount(); i++)
            {
                for (size_t j = 0; j < m_Layers[i].numNodes; j++)
                {
                    Nodes.emplace_back(m_Layers[i].material);
                }
                NumNodes.push_back(m_Layers[i].numNodes);
            }

            std::partial_sum(NumNodes.begin(), NumNodes.end(), NumNodes.begin());
            NodesAreValid = true;
        }

        // For initializing temperatures to layer stack
        friend class ResponseSolverBase;
    };

} // namespace OTAP
