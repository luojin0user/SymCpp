#pragma once
#include <cassert>
#include <algorithm>
#include <cstring>

#include <Eigen/Dense>

class Matrix2D
{
public:
    Matrix2D() = default;

    Matrix2D(size_t rows, size_t cols)
        : rows_(rows), cols_(cols)
    {
        data_ = new double[rows_ * cols_]();
    }

    ~Matrix2D()
    {
        delete[] data_;
    }

    // 禁止拷贝（避免误拷贝大矩阵）
    Matrix2D(const Matrix2D &) = delete;
    Matrix2D &operator=(const Matrix2D &) = delete;

    // 允许移动（所有权转移）
    Matrix2D(Matrix2D &&other) noexcept
    {
        move_from(other);
    }

    Matrix2D &operator=(Matrix2D &&other) noexcept
    {
        if (this != &other)
        {
            delete[] data_;
            move_from(other);
        }
        return *this;
    }

    Matrix2D &operator+=(const Matrix2D &rhs)
    {
        assert(rows_ == rhs.rows_);
        assert(cols_ == rhs.cols_);

        size_t n = rows_ * cols_;

        // OpenMP 并行（如果你已经开启）
        // #pragma omp parallel for
        for (size_t i = 0; i < n; ++i)
        {
            data_[i] += rhs.data_[i];
        }

        return *this;
    }

    /* ===== 单元素访问 ===== */
    inline double &at(size_t i, size_t j)
    {
        assert(i < rows_ && j < cols_);
        return data_[i * cols_ + j];
    }

    inline const double &at(size_t i, size_t j) const
    {
        assert(i < rows_ && j < cols_);
        return data_[i * cols_ + j];
    }

    /* ===== 子块写入 =====
       src:     外部连续存储的数据（行主序）
       src_rows × src_cols
       写入到 [sr, er), [sc, ec)
    */
    void write_block(
        size_t start_row,
        size_t start_col,
        const Matrix2D &block)
    {
        // 边界检查
        assert(start_row + block.rows_ <= rows_);
        assert(start_col + block.cols_ <= cols_);

        for (size_t r = 0; r < block.rows_; ++r)
        {
            std::memcpy(
                data_ + (start_row + r) * cols_ + start_col,
                block.data_ + r * block.cols_,
                block.cols_ * sizeof(double));
        }
    }

    /* ===== 子块读取 ===== */
    void read_block(
        size_t sr, size_t sc,
        size_t er, size_t ec,
        Matrix2D &dst) const
    {
        assert(sr < er && sc < ec);
        assert(er <= rows_ && ec <= cols_);

        size_t block_rows = er - sr;
        size_t block_cols = ec - sc;

        // 目标矩阵尺寸必须匹配子块
        assert(dst.rows() == block_rows);
        assert(dst.cols() == block_cols);

        for (size_t r = 0; r < block_rows; ++r)
        {
            std::memcpy(
                dst.data() + r * block_cols,
                data_ + (sr + r) * cols_ + sc,
                block_cols * sizeof(double));
        }
    }

    /* ===== 填充 ===== */
    void fill(double value)
    {
        std::fill(data_, data_ + rows_ * cols_, value);
    }

    Eigen::MatrixXd toEigenCopy() const
    {
        return Eigen::Map<
            const Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>>(data_, rows_, cols_);
    }

    /* ===== 基本信息 ===== */
    size_t rows() const { return rows_; }
    size_t cols() const { return cols_; }
    double *data() { return data_; }
    const double *data() const { return data_; }

private:
    double *data_ = nullptr;
    size_t rows_ = 0;
    size_t cols_ = 0;

    void move_from(Matrix2D &other)
    {
        data_ = other.data_;
        rows_ = other.rows_;
        cols_ = other.cols_;
        other.data_ = nullptr;
        other.rows_ = other.cols_ = 0;
    }
};
