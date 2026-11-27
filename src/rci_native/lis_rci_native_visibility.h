/* Internal ELF visibility; qualification builds may explicitly expose hooks. */
#ifndef LIS_RCI_NATIVE_VISIBILITY_H
#define LIS_RCI_NATIVE_VISIBILITY_H
#if defined(__GNUC__) && !defined(_WIN32) && !defined(LIS_RCI_NATIVE_TEST_HOOKS)
#define LIS_RCI_PRIVATE __attribute__((visibility("hidden")))
#else
#define LIS_RCI_PRIVATE
#endif
#endif
