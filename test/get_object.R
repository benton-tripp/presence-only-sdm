# Get model or other object from cache if it has been saved before
get.object <- function(obj, file.name, obj.path, read.func=readRDS, save.func=saveRDS, ...) {
  f.path <- file.path(obj.path, file.name)
  if (!dir.exists(obj.path)) {
    dir.create(obj.path)
  }
  # Object cache check
  if (file.exists(f.path)) {
    obj <- read.func(f.path)
  } else {
    save.func(obj, f.path, ...)
  }
  obj
}