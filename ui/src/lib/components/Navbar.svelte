<script lang="ts">
 import * as NavigationMenu from "$lib/components/ui/navigation-menu/index.js";
</script>

<div class="px-4 py-5 content-center flex justify-center">
<NavigationMenu.Root>
 <NavigationMenu.List>
  <NavigationMenu.Item>
   <NavigationMenu.Trigger>Item One</NavigationMenu.Trigger>
   <NavigationMenu.Content>
    <NavigationMenu.Link>Link</NavigationMenu.Link>
   </NavigationMenu.Content>
  </NavigationMenu.Item>
 </NavigationMenu.List>
</NavigationMenu.Root>
</div>
