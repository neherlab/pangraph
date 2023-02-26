export type SetterOrUpdater<T> = (valOrUpdater: ((currVal: T) => T) | T) => void
