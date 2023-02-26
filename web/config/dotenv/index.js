import { join } from 'path'
import { config } from 'dotenv'
import { findModuleRoot } from '../../lib/findModuleRoot'

const { moduleRoot } = findModuleRoot()
config({ path: join(moduleRoot, '.env') })
