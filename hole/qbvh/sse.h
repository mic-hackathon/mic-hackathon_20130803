#ifndef _SSE_H_
#define _SSE_H_

#ifdef _MSC_VER
    #include <malloc.h>
#else
    #include <stdlib.h>
    static inline void *_aligned_malloc(size_t size, size_t alignment)
    {
        void *p;
        int ret = posix_memalign(&p, alignment, size);
        return (ret == 0) ? p : 0;
    }
#endif
    /* 16byteアライメントされた500byteのデータを確保 */
    // char *p = (char*)_aligned_malloc(500, 16);

#endif
