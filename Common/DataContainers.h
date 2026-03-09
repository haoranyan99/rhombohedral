#pragma once
#include <vector>
#include <string>
#include <memory>
#include <stdexcept>
#include <Eigen/Dense>
#include <unordered_map>   // NEW
#include <cstdint>         // NEW (for int64_t)


namespace core {

// ============================================================
// Field system (heterogeneous typed fields)
// ============================================================
struct FieldBase {
    std::string name;
    virtual ~FieldBase() = default;
    virtual size_t size() const = 0;
    virtual void resize(size_t N) = 0;   // NEW
};

template <class T>
struct Field : FieldBase {
    std::vector<T> v;

    Field(const std::string& n, size_t N) {
        name = n;
        v.resize(N);
    }

    size_t size() const override { return v.size(); }
    void resize(size_t N) override { v.resize(N); }   // NEW
};


// ============================================================
// Helper: add/get typed fields (shared by GridData / PathDataT)
// ============================================================
namespace detail {

inline void ensure_field_not_exists_(const std::vector<std::shared_ptr<FieldBase>>& fields,
                                     const std::string& name,
                                     const char* where)
{
    for (auto& f : fields) {
        if (f->name == name) {
            throw std::runtime_error(std::string(where) + ": field already exists: " + name);
        }
    }
}

template <class T>
inline Field<T>& get_field_(std::vector<std::shared_ptr<FieldBase>>& fields,
                           const std::string& name,
                           const char* where)
{
    for (auto& f : fields) {
        if (f->name == name) {
            auto* p = dynamic_cast<Field<T>*>(f.get());
            if (!p) {
                throw std::runtime_error(std::string(where) + ": type mismatch for field: " + name);
            }
            return *p;
        }
    }
    throw std::runtime_error(std::string(where) + ": field not found: " + name);
}

template <class T>
inline const Field<T>& get_field_(const std::vector<std::shared_ptr<FieldBase>>& fields,
                                 const std::string& name,
                                 const char* where)
{
    for (auto& f : fields) {
        if (f->name == name) {
            auto* p = dynamic_cast<const Field<T>*>(f.get());
            if (!p) {
                throw std::runtime_error(std::string(where) + ": type mismatch for field: " + name);
            }
            return *p;
        }
    }
    throw std::runtime_error(std::string(where) + ": field not found: " + name);
}

} // namespace detail







// ============================================================
// GridData: 2D integer grid container
// ------------------------------------------------------------
// - iq/jq: integer coordinates (flattened list)
// - dx/dy: generic step parameters (interpretation is context-dependent)
// - fields: arbitrary per-grid-point data
// ============================================================
struct GridData {
    std::vector<int> iq;
    std::vector<int> jq;

    double dx = 0.0;
    double dy = 0.0;

    std::vector<std::shared_ptr<FieldBase>> fields;

    // --------------------------
    // NEW: (iq,jq) -> idx cache
    // --------------------------
    mutable bool idx_cache_valid = false;
    mutable std::unordered_map<long long, size_t> ij2idx;

    size_t size() const { return iq.size(); }

    void resize(size_t N) {
        iq.resize(N);
        jq.resize(N);
        for (auto& f : fields) f->resize(N);

        // NEW: invalidate cache on resize
        idx_cache_valid = false;
        ij2idx.clear();
    }

    void assert_consistent() const {
        const size_t N = size();
        if (jq.size() != N)
            throw std::runtime_error("GridData: iq/jq size mismatch");
        for (auto& f : fields) {
            if (f->size() != N)
                throw std::runtime_error("GridData: field size mismatch: " + f->name);
        }
    }

    template <class T>
    Field<T>& add(const std::string& name) {
        detail::ensure_field_not_exists_(fields, name, "GridData::add");
        auto p = std::make_shared<Field<T>>(name, size());
        fields.push_back(p);
        return *p;
    }

    template <class T>
    Field<T>& get(const std::string& name) {
        return detail::get_field_<T>(fields, name, "GridData::get");
    }

    template <class T>
    const Field<T>& get(const std::string& name) const {
        return detail::get_field_<T>(fields, name, "GridData::get");
    }

    // ============================================================
    // NEW helper 1: idx -> (iq,jq)  (no dependence on fill order)
    // ============================================================
    inline std::pair<int,int> idx_to_ij(size_t idx) const {
        if (idx >= size())
            throw std::runtime_error("GridData::idx_to_ij: idx out of range");
        return { iq[idx], jq[idx] };
    }

    // ============================================================
    // NEW helper 2: (iq,jq) -> idx, via lazy-built cache
    // - Returns true if found; false if this (iq,jq) doesn't exist in grid
    // ============================================================
    inline bool ij_to_idx(int i, int j, size_t& out_idx) const {
        build_ij2idx_cache_if_needed_();
        const long long key = pack_ij_(i, j);
        auto it = ij2idx.find(key);
        if (it == ij2idx.end()) return false;
        out_idx = it->second;
        return true;
    }

    // convenience: throw-on-miss version
    inline size_t ij_to_idx_or_throw(int i, int j) const {
        size_t idx;
        if (!ij_to_idx(i, j, idx)) {
            throw std::runtime_error("GridData::ij_to_idx: (iq,jq) not found in grid");
        }
        return idx;
    }

private:
    // pack two 32-bit ints into one 64-bit key (stable, fast, no struct)
    static inline long long pack_ij_(int i, int j) {
        return ( (long long)( (uint32_t)i ) << 32 ) | (uint32_t)j;
    }

    inline void build_ij2idx_cache_if_needed_() const {
        if (idx_cache_valid) return;

        ij2idx.clear();
        ij2idx.reserve(size() * 2);

        for (size_t idx = 0; idx < size(); ++idx) {
            ij2idx[ pack_ij_(iq[idx], jq[idx]) ] = idx;
        }

        idx_cache_valid = true;
    }
};






// ============================================================
// PathData: 1D k-path container (non-template)
// ------------------------------------------------------------
// - k_list: k points along the path (Eigen::Vector2d)
// - kline : cumulative distance along the path
// - xticks: tick positions and labels for plotting
// - fields: arbitrary per-k-point data
// ============================================================
struct PathData {
    using Vec2 = Eigen::Vector2d;

    std::vector<Vec2>   k_list;
    std::vector<double> kline;

    std::vector<double>      xtick_pos;
    std::vector<std::string> xtick_lab;

    std::vector<std::shared_ptr<FieldBase>> fields;

    size_t size() const { return k_list.size(); }

    void resize(size_t N) {
        k_list.resize(N);
        kline.resize(N);
        for (auto& f : fields) f->resize(N);
    }

    void assert_consistent() const {
        const size_t N = size();
        for (auto& f : fields) {
            if (f->size() != N)
                throw std::runtime_error("PathData: field size mismatch: " + f->name);
        }
    }

    template <class T>
    Field<T>& add(const std::string& name) {
        detail::ensure_field_not_exists_(fields, name, "PathData::add");
        auto p = std::make_shared<Field<T>>(name, size());
        fields.push_back(p);
        return *p;
    }

    template <class T>
    Field<T>& get(const std::string& name) {
        return detail::get_field_<T>(fields, name, "PathData::get");
    }

    template <class T>
    const Field<T>& get(const std::string& name) const {
        return detail::get_field_<T>(fields, name, "PathData::get");
    }
};







// ============================================================
// SeriesData: 1D x-axis container (generic 1D series)
// ------------------------------------------------------------
// - x: 1D axis (e.g., energy, frequency, filling, etc.)
// - dx: optional uniform step (0 means "unknown / non-uniform")
// - fields: arbitrary per-point data (same length as x)
// ============================================================
struct SeriesData {
    std::vector<double> x;   // axis
    double dx = 0.0;         // optional step; can be 0 for non-uniform

    std::vector<std::shared_ptr<FieldBase>> fields;

    size_t size() const { return x.size(); }

    void resize(size_t N) {
        x.resize(N);
        for (auto& f : fields) f->resize(N);
    }

    void assert_consistent() const {
        const size_t N = size();
        for (auto& f : fields) {
            if (f->size() != N)
                throw std::runtime_error("SeriesData: field size mismatch: " + f->name);
        }
    }

    template <class T>
    Field<T>& add(const std::string& name) {
        detail::ensure_field_not_exists_(fields, name, "SeriesData::add");
        auto p = std::make_shared<Field<T>>(name, size());
        fields.push_back(p);
        return *p;
    }

    template <class T>
    Field<T>& get(const std::string& name) {
        return detail::get_field_<T>(fields, name, "SeriesData::get");
    }

    template <class T>
    const Field<T>& get(const std::string& name) const {
        return detail::get_field_<T>(fields, name, "SeriesData::get");
    }
};




























} // namespace core
