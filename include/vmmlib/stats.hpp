/* 
 * File:   result.hpp
 * Author: rballester
 *
 * Created on August 29, 2012, 3:13 PM
 */
#include <string>

#ifndef RESULT_HPP
#define	RESULT_HPP

namespace vmml {

    class stats {
    public:
        stats();
        stats(const stats& other);
        inline stats operator+(const stats& other) const;
        void operator+=(const stats& other);

        friend std::ostream& operator <<(std::ostream& os,
                const stats& stats) {
            os << stats.description << " " << stats.ranks << " " << stats.n_iterations << " " << stats.error << " " << stats.nnz << " " << stats.dec_time << " " << stats.rec_time;
            return os;
        }

        std::string get_description();
        std::string get_short_description();
        int get_ranks();
        int get_n_iterations();
        double get_error();
        int get_nnz();
        double get_dec_time();
        double get_rec_time();

        void set_description(std::string description);
        void set_short_description(std::string short_description);
        void set_ranks(int ranks);
        void set_n_iterations(int n_iterations);
        void set_error(double error);
        void set_nnz(int nnz);
        void set_dec_time(double dec_time);
        void set_rec_time(double rec_time);

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

    stats::stats() {
        ranks = n_iterations = error = nnz = dec_time = rec_time = 0;
    }

    stats::stats(const stats& other) {
        description = other.description;
        ranks = other.ranks;
        n_iterations = other.n_iterations;
        error = other.error;
        nnz = other.nnz;
        dec_time = other.dec_time;
        rec_time = other.rec_time;
    }

    inline stats stats::operator+(const stats& other) const {
        stats result(*this);
        result += other;
        return result;
    }

    void stats::operator+=(const stats& other) {
        ranks += other.ranks;
        n_iterations += other.n_iterations;
        error += other.error;
        nnz += other.nnz;
        dec_time += other.dec_time;
        rec_time += other.rec_time;
    }

    std::string stats::get_description() {
        return description;
    }
    
    std::string stats::get_short_description() {
        return short_description;
    }

    int stats::get_ranks() {
        return ranks;
    }
    
    int stats::get_n_iterations() {
        return n_iterations;
    }

    double stats::get_error() {
        return error;
    }

    int stats::get_nnz() {
        return nnz;
    }

    double stats::get_dec_time() {
        return dec_time;
    }

    double stats::get_rec_time() {
        return rec_time;
    }

    void stats::set_description(std::string description) {
        this->description = description;
    }
    
    void stats::set_short_description(std::string short_description) {
        this->short_description = short_description;
    }

    void stats::set_ranks(int ranks) {
        this->ranks = ranks;
    }
    
    void stats::set_n_iterations(int n_iterations) {
        this->n_iterations = n_iterations;
    }

    void stats::set_error(double error) {
        this->error = error;
    }

    void stats::set_nnz(int nnz) {
        this->nnz = nnz;
    }

    void stats::set_dec_time(double dec_time) {
        this->dec_time = dec_time;
    }

    void stats::set_rec_time(double rec_time) {
        this->rec_time = rec_time;
    }

} // namespace vmml

#endif	/* RESULT_HPP */

