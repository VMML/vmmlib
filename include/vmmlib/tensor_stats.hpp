/* 
 * VMMLib - Tensor Classes
 *  
 * @author Rafael Ballester
 *
 * Class for encapsulating several tensor experimental properties, so that they can be easily returned, manipulated and printed (approximation error, decomposition and reconstruction times, etc.).
 * 
 */
#include <string>

#ifndef __VMML__TENSOR_STATS__HPP__
#define	__VMML__TENSOR_STATS__HPP__

namespace vmml {

    class tensor_stats {
    public:

        tensor_stats() {
            ranks = n_iterations = error = nnz = dec_time = rec_time = 0;
        }

        tensor_stats(const tensor_stats& other) {
            description = other.description;
            ranks = other.ranks;
            n_iterations = other.n_iterations;
            error = other.error;
            nnz = other.nnz;
            dec_time = other.dec_time;
            rec_time = other.rec_time;
        }

        inline tensor_stats operator+(const tensor_stats& other) const {
            tensor_stats result(*this);
            result += other;
            return result;
        }

        void operator+=(const tensor_stats& other) {
            ranks += other.ranks;
            n_iterations += other.n_iterations;
            error += other.error;
            nnz += other.nnz;
            dec_time += other.dec_time;
            rec_time += other.rec_time;
        }

        friend std::ostream& operator <<(std::ostream& os,
                const tensor_stats& stats) {
            os << stats.description << " " << stats.ranks << " " << stats.n_iterations << " " << stats.error << " " << stats.nnz << " " << stats.dec_time << " " << stats.rec_time;
            return os;
        }

        std::string get_description() {
            return description;
        }

        std::string get_short_description() {
            return short_description;
        }

        int get_ranks() {
            return ranks;
        }

        int get_n_iterations() {
            return n_iterations;
        }

        double get_error() {
            return error;
        }

        int get_nnz() {
            return nnz;
        }

        double get_dec_time() {
            return dec_time;
        }

        double get_rec_time() {
            return rec_time;
        }

        void set_description(std::string description) {
            this->description = description;
        }

        void set_short_description(std::string short_description) {
            this->short_description = short_description;
        }

        void set_ranks(int ranks) {
            this->ranks = ranks;
        }

        void set_n_iterations(int n_iterations) {
            this->n_iterations = n_iterations;
        }

        void set_error(double error) {
            this->error = error;
        }

        void set_nnz(int nnz) {
            this->nnz = nnz;
        }

        void set_dec_time(double dec_time) {
            this->dec_time = dec_time;
        }

        void set_rec_time(double rec_time) {
            this->rec_time = rec_time;
        }

        static std::string get_content() {
            return "description n_iterations error nnz dec_time rec_time";
        }

    private:
        std::string description;
        std::string short_description; // The description, without rank sizes
        int ranks;
        int n_iterations;
        double error;
        int nnz;
        double dec_time;
        double rec_time;

    };

} // namespace vmml

#endif

