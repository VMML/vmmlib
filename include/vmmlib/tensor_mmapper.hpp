#ifndef __VMML__TENSOR_MMAPPER__HPP__
#define __VMML__TENSOR_MMAPPER__HPP__

#include "tensor4.hpp"
#include <sys/mman.h>


namespace vmml {

    template< typename T > // ein tensor3 oder tensor4 typ
    class tensor_mmapper {
    public:
        //tensor_mapper( const std::string& dir_, const std::string& filename, size_t tensor_count = 1 );

        //allocate memory mapped file
        tensor_mmapper(const std::string& dir_, const std::string& filename_, const bool prot_read_);
        ~tensor_mmapper();

        void get_tensor(T& tensor_) {
            tensor_ = *_tensor;
        }

    protected:
        static std::string concat_path(const std::string& dir_, const std::string& filename_);

        void* t3_allocate_rd_mmap(const std::string& dir_, const std::string& filename_);
        void* t3_allocate_wr_mmap(const std::string& dir_, const std::string& filename_);

        int _fd; // one mmap for all tensors
        size_t _file_size; // sizeof tensors * tensor_count
        size_t _tensor_size;
        void* _data; // ptr to mmapped memory 

        T* _tensor;

        //std::vector<T*> _tensors;

    }; //end tensor_mmapper



#define VMML_TEMPLATE_STRING        template< typename T >
#define VMML_TEMPLATE_CLASSNAME     tensor_mmapper< T >

    /*VMML_TEMPLATE_STRING
    VMML_TEMPLATE_CLASSNAME::tensor_mmapper( const std::string& dir_, const std::string& filename, size_t tensor_count = 1 ) 
    // : init all values 
    {
            _file_size = ... 
            // open or create file and mmap it
	
            _data = ...
            _tensor_size = ... 
	
            void* d = _data;
            for( size_t i = 0; i < tensor_count; ++i ) 
            {
                    _tensors.push_back( new T(d));
                    d += _tensor_size;
            }
    }*/


    VMML_TEMPLATE_STRING
    VMML_TEMPLATE_CLASSNAME::~tensor_mmapper() {
        if (_tensor) {
            _tensor->clear_array_pointer();
            delete _tensor;
        }

        munmap(_data, _file_size); //get error
        close(_fd);

    }

    VMML_TEMPLATE_STRING
    VMML_TEMPLATE_CLASSNAME::tensor_mmapper(const std::string& dir_, const std::string& filename_, const bool prot_read_)
    : _file_size(0), _data(0) {

        _file_size = T::get_array_size_in_bytes();

        if (prot_read_) {
            t3_allocate_rd_mmap(dir_, filename_);
        } else {
            t3_allocate_wr_mmap(dir_, filename_);
        }

        if (_data != 0) {
            _tensor = new T(_data);
        }
    }

    VMML_TEMPLATE_STRING
    void*
    VMML_TEMPLATE_CLASSNAME::
    t3_allocate_rd_mmap(const std::string& dir_, const std::string& filename_) {
        //void * mmap(void *addr, size_t len, int prot, int flags, int fd, off_t offset);

        std::string path = concat_path(dir_, filename_);

        _fd = open(path.c_str(), O_RDONLY, 0);
        if (_fd == -1) {
            {
                close(_fd);
                std::cout << "no file open for memory mapping" << std::endl;
            }
        }


        off_t offset = 0;

        _data = (void*) mmap(0, _file_size, PROT_READ, MAP_FILE | MAP_SHARED, _fd, offset);

        if (_data == MAP_FAILED) {
            std::cout << "mmap failed" << std::endl;
            return 0;
        }

        return _data;
    }

    VMML_TEMPLATE_STRING
    void*
    VMML_TEMPLATE_CLASSNAME::
    t3_allocate_wr_mmap(const std::string& dir_, const std::string& filename_) {
        //void * mmap(void *addr, size_t len, int prot, int flags, int fd, off_t offset);

        std::string path = concat_path(dir_, filename_);

        _fd = open(path.c_str(), O_RDWR, 0);
        if (_fd == -1) {
            if (!_tensor) {
                _tensor = new T();
            }
            _tensor->zero();
            //                TODO Susanne: can you translate this line?
            //		_tensor->write_to_raw( dir_, filename_ );

            _fd = open(path.c_str(), O_RDWR, 0);
            if (_fd == -1) {
                close(_fd);
                std::cout << "no file open for memory mapping" << std::endl;
            }
        }


        off_t offset = 0;

        _data = (void*) mmap(0, _file_size, PROT_WRITE, MAP_FILE | MAP_SHARED, _fd, offset);

        if (_data == MAP_FAILED) {
            std::cout << "mmap failed" << std::endl;
            return 0;
        }

        return _data;
    }

    VMML_TEMPLATE_STRING
    std::string
    VMML_TEMPLATE_CLASSNAME::concat_path(const std::string& dir_, const std::string& filename_) {
        std::string path = "";
        int dir_length = dir_.size() - 1;
        int last_separator = dir_.find_last_of("/");
        path = dir_;
        if (last_separator < dir_length) {
            path.append("/");
        }
        path.append(filename_);

        //check for format
        if (filename_.find("raw", filename_.size() - 3) == (-1)) {
            path.append(".");
            path.append("raw");
        }

        return path;
    }



#undef VMML_TEMPLATE_STRING
#undef VMML_TEMPLATE_CLASSNAME


}//end vmml namespace

#endif
